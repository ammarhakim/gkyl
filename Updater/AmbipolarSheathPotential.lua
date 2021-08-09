-- Gkyl ------------------------------------------------------------------------
--
-- Asumming ambipolar fluxes at the sheath entrance, compute the sheath
-- and the electrostatic potential in the whole domain.
-- The potential is computed from
--    phi = phi_s - (T_e/q_e)*log(n_e - n_{es})
-- where phi_s and n_{es} are the sheath potential and electron density,
-- respectively. Since we typically have two sheaths, we compute this in two
-- parts. Continuity is later enforced by an FEM smoothing operator. The sheath
-- potential is
--    phi_s = (T_e/q_e)*log( sqrt(2*pi)*Gamma_i/(n_e*v_{te}) )
-- with Gamma_i being the ion particle flux through the sheath entrance.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase  = require "Updater.Base"
local DataStruct   = require "DataStruct"
local LinearDecomp = require "Lib.LinearDecomp"
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local xsys         = require "xsys"
local ModDecl      = require "Updater.ambipolarSheathPotentialData.AmbipolarSheathPotentialModDecl"
local Updater      = require "Updater"

-- Inherit the base Updater from UpdaterBase updater object.
local ASPotential = Proto(UpdaterBase)

function ASheath:init(tbl)
   ASheath.super.init(self, tbl) -- setup base object

   self.grid  = assert(tbl.onGrid, "Updater.AmbipolarSheathPotential: Must provide configuration space grid object 'onGrid'.")
   self.basis = assert(tbl.basis, "Updater.AmbipolarSheathPotential: Must provide configuration space basis object 'basis'.")
   self.boundaryGrids = assert(tbl.basis, "Updater.AmbipolarSheathPotential: Must provide configuration space boundary grids 'boundaryGrids'.")

   self.mElc  = assert(tbl.electronMass, "Updater.AmbipolarSheathPotential: Must provide electron mass in 'electronMass'.")
   self.qElc  = assert(tbl.electronCharge, "Updater.AmbipolarSheathPotential: Must provide electron charge in 'electronCharge'.")

   -- Number of quadrature points in each direction
   local numQuad1D = self.basis:polyOrder() + 1

   self.quadType = "Gauss"

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self.dim = self.basis:ndim()
   self.sheathDir = self.dim   -- Assume the sheath direction is the last dimension.

   self._sheathEvKer = ModDecl.selectPhiSheathQuad(self.basis:id(), self.dim, self.basis:polyOrder(), self.quadType)
   self._phiKer      = ModDecl.selectPhiQuad(self.basis:id(), self.dim, self.basis:polyOrder(), self.quadType)

   self.idxB = Lin.IntVec(self.dim)
   local xcB = Lin.Vec(self.dim)

   -- We compute the potential in two halfs [z_min,0] and [0,z_max] first.
   -- For this we need two pairs of grids.
   --   - A lower/upper pair of skin-cell grids to compute phi_s (boundaryGrids).
   --   - The [z_min,0] and [0,z_max] to compute phi (halfGrid).
   self.boundary = {"lower","upper"}
   self.halfSign = {lower=-1., upper=1.}
   self.halfDomRangeDecomp = {}
   local sheathDirCells = self.grid:numCells(self.sheathDir)
   local globalRange = self.grid:globalRange()
   for _, b in ipairs(self.boundary) do
      -- Compute the offsets used to shorten the domain in the sheath direction. Assume that the
      -- zero in the sheath direction is located at a cell boundary and not inside of a cell.
      local sheathDirExts = {0,0}
      for idx = 1, sheathDirCells do
         self.idxB[self.sheathDir] = idx
         self.grid:setIndex(self.idxB)
         self.grid:cellCenter(xcB)
         if (b == "upper") and (xcB[self.sheathDir] > 0.0) then
            sheathDirExts[1] = -(idx-1)
            break
         elseif (b == "lower") and (xcB[self.sheathDir] > 0.0) then
            sheathDirExts[2] = -(sheathDirCells-(idx-1))
            break
         end
      end
      -- Create the decomposed range that a rank has to
      -- loop over when looping over the half domain.
      local halfDomRange = self.grid:globalRange():extendDir(self.sheathDir,sheathDirExts[1],sheathDirExts[2])
      local halfDomLocalRange = self.grid:localRange():intersect(halfDomRange)
      self.halfDomRangeDecomp[b] = LinearDecomp.LinearDecompRange {
        range = halfDomLocalRange, numSplit = self.grid:numSharedProcs() }
   end

   self.phiSheath = {lower=DataStruct.Field{onGrid = self.boundaryGrids["lower"],
                                            numComponents = self.basis:numBasis(),
                                            ghost = {1,1},},
                     upper=DataStruct.Field{onGrid = self.boundaryGrids["upper"], 
                                            numComponents = self.basis:numBasis(),
                                            ghost = {1,1},}}
   self.m0ElcSheath = {lower=DataStruct.Field{onGrid = self.boundaryGrids["lower"],
                                              numComponents = self.basis:numBasis(),
                                              ghost = {1,1},},
                       upper=DataStruct.Field{onGrid = self.boundaryGrids["upper"], 
                                              numComponents = self.basis:numBasis(),
                                              ghost = {1,1},}}
   -- Pre-define some pointers and indexers.
   self.m0ElcSheathPtr = {lower=0, upper=0}
   self.phiSheathPtr   = {lower=0, upper=0}
   self.GammaIonPtr    = {lower=0, upper=0}   -- Set in :advance.
   self.boundaryIdxr   = {lower=0, upper=0}
   for _, b in ipairs(self.boundary) do
      self.phiSheath[b]:clear(0.0)
      self.m0ElcSheath[b]:clear(0.0)

      self.phiSheathPtr[b]   = phiSheath[b]:get(1)
      self.m0ElcSheathPtr[b] = m0ElcSheath[b]:get(1)

      self.boundaryIdxr[b] = self.phiSheath[b]:genIndexer()
   end

   self.phiZSmoother = Updater.FemParPoisson {
      onGrid = self.grid,   bcLower = {{T="N",V=0.0}},
      basis  = self.basis,  bcUpper = {{T="N",V=0.0}},
      smooth = true,
   }
   self.phiAux = DataStruct.Field{onGrid = self.grid,
                                  numComponents = self.basis:numBasis(),
                                  ghost = {1,1},}
   self.phiAux:clear(0.)
   self.phiAuxPtr = self.phiAux:get(1)
end

function ASheath:_advance(tCurr, inFlds, outFlds)
   local GammaIon, m0Elc, vtSqElc = inFlds[1], inFlds[2], inFlds[3]
   local phi = outFlds[1]

   self.GammaIonPtr["lower"]  = GammaIon["lower"]:get(1)
   self.GammaIonPtr["upper"]  = GammaIon["upper"]:get(1)
   local m0ElcPtr, vtSqElcPtr = m0Elc:get(1), vtSqElc:get(1)
   local phiPtr               = phi:get(1)
   local domainIdxr           = inFld:genIndexer()

   local grid = self.grid
   local globalRange = phi:globalRange()

   for _, b in ipairs(self.boundary) do
      -- Loop over boundary grid and compute phi_s in each cell using quadrature.
      local skinRange = b=="lower" and globalRange:lowerSkin(self.sheathDir,1) or globalRange:upperSkin(self.sheathDir,1)
      local skinRangeDecomp = LinearDecomp.LinearDecompRange {
            range = skinRange, numSplit = grid:numSharedProcs() }
      local tId = grid:subGridSharedId()    -- Local thread ID.
      for idx in skinRangeDecomp:rowMajorIter(tId) do
         idx:copyInto(self.idxB)
         self.idxB[self.sheathDir] = 1

         GammaIon[b]:fill(boundaryIdxr[b](self.idxB), GammaIonPtr[b])
         phiSheath:fill(boundaryIdxr(self.idxB), phiSheathPtr)
         m0ElcSheath:fill(boundaryIdxr(self.idxB), m0ElcSheathPtr)
         m0Elc:fill(innerIdxr(idx), m0ElcPtr)
         vtSqElc:fill(innerIdxr(idx), vtSqElcPtr)
         
         self._sheathEvKer[b](qElc, mElc, GammaIonPtr[b]:data(), m0ElcPtr:data(), vtSqElcPtr:data(), m0ElcSheathPtr:data(), phiSheathPtr:data())
      end

      -- Loop over the half grid and compute phi.
      for idx in self.halfDomRangeDecomp[b]:rowMajorIter(tId) do
         idx:copyInto(self.idxB)
         self.idxB[self.sheathDir] = 1

         phiSheath:fill(boundaryIdxr(self.idxB), phiSheathPtr)
         m0ElcSheath:fill(boundaryIdxr(self.idxB), m0ElcSheathPtr)
         vtSqElc:fill(innerIdxr(idx), vtSqElcPtr)
         m0Elc:fill(innerIdxr(idx), m0ElcPtr)
         self.phiAux:fill(innerIdxr(idx), self.phiAuxPtr)

         self._phiKer(self.halfSign[b], qElc, mElc, m0ElcPtr:data(), vtSqElcPtr:data(), m0ElcSheathPtr:data(), phiSheathPtr:data(), self.phiAuxPtr:data())
      end
   end

   -- Use the FEM smoothing operator to enforce continuity.
   self.phiZSmoother:advance(tCurr, {self.phiAux}, {phi})
end

return ASheath
