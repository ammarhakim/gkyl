-- Gkyl ------------------------------------------------------------------------
--
-- Ambipolar, adiabatic electron Sheath Potential (ASheathPotential):
-- Asumming ambipolar fluxes at the sheath entrance, compute the sheath
-- and the electrostatic potential in the whole domain.
-- The potential is computed from
--    phi = phi_s - (T_e/q_e)*log(n_e/n_{es})
-- which thanks to quasineutrality becomes
--    phi = phi_s - (T_e/q_e)*log(n_i/n_{is})
-- where phi_s and n_{ks} are the potential and density of species k at the
-- sheath entrance, respectively. Since we typically have two sheaths, we
-- compute this in two parts. Continuity is later enforced by an FEM smoothing
-- operator. The sheath potential is
--    phi_s = -(T_e/q_e)*log( sqrt(2*pi)*Gamma_i/(n_e*v_{te}) )
-- which using quasineutrality again gives
--    phi_s = -(T_e/q_e)*( log( sqrt(2*pi)*Gamma_i/(n_i*v_{te}) ) )|_{z=z_s}
-- with Gamma_i being the ion particle flux, and the natural logarithm being
-- evaluated at the sheath entrance.
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
local ModDecl      = require "Updater.aSheathPotentialData.asheath_potential_mod_decl"
local Mpi          = require "Comm.Mpi"

-- Inherit the base Updater from UpdaterBase updater object.
local ASheathPotential = Proto(UpdaterBase)

function ASheathPotential:init(tbl)
   ASheathPotential.super.init(self, tbl) -- setup base object

   self.grid  = assert(tbl.onGrid, "Updater.ASheathPotential: Must provide configuration space grid object 'onGrid'.")
   self.basis = assert(tbl.basis, "Updater.ASheathPotential: Must provide configuration space basis object 'basis'.")
   self.boundaryGrids = assert(tbl.boundaryGrids, "Updater.ASheathPotential: Must provide configuration space boundary grids 'boundaryGrids'.")

   self.mElc    = assert(tbl.electronMass, "Updater.ASheathPotential: Must provide electron mass in 'electronMass'.")
   self.qElc    = assert(tbl.electronCharge, "Updater.ASheathPotential: Must provide electron charge in 'electronCharge'.")
   self.tempElc = assert(tbl.electronTemp, "Updater.ASheathPotential: Must provide electron temperature in 'electronTemp'.")

   -- Number of quadrature points in each direction
   local numQuad1D = self.basis:polyOrder() + 1

   self.quadType = "Gauss"

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self.dim = self.basis:ndim()
   self.sheathDir = self.dim   -- Assume the sheath direction is the last dimension.

   self._phiSheathKer = ModDecl.selectPhiSheathQuad(self.basis:id(), self.dim, self.basis:polyOrder(), self.quadType)
   self._phiKer       = ModDecl.selectPhiQuad(self.basis:id(), self.dim, self.basis:polyOrder(), self.quadType)

   self.dx   = Lin.Vec(self.dim)
   self.idxB = Lin.IntVec(self.dim)
   for d=1,self.dim do self.idxB[d] = 1 end
   local xcB = Lin.Vec(self.dim)

   -- We compute the potential in two halfs [z_min,0] and [0,z_max] first.
   -- For this we need two pairs of grids.
   --   - A lower/upper pair of skin-cell grids to compute phi_s (boundaryGrids).
   --   - The [z_min,0] and [0,z_max] grids to compute phi (halfGrid).
   self.boundary = {"lower","upper"}
   self.halfSign = {lower=-1., upper=1.}
   self.halfDomRangeDecomp = {}
   local sheathDirCells = self.grid:numCells(self.sheathDir)
   local globalRange = self.grid:globalRange()
   for _, b in ipairs(self.boundary) do
      -- Compute the offsets used to shorten the domain in the sheath direction.
      -- ASSUME zero in the sheath direction is at a cell boundary and not inside of a cell.
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
   self.m0IonSheath = {lower=DataStruct.Field{onGrid = self.boundaryGrids["lower"],
                                              numComponents = self.basis:numBasis(),
                                              ghost = {1,1},},
                       upper=DataStruct.Field{onGrid = self.boundaryGrids["upper"], 
                                              numComponents = self.basis:numBasis(),
                                              ghost = {1,1},}}
   -- Pre-define some pointers and indexers.
   self.m0IonSheathPtr = {lower=0, upper=0}
   self.phiSheathPtr   = {lower=0, upper=0}
   self.GammaIonPtr    = {lower=0, upper=0}   -- Set in :advance.
   self.boundaryIdxr   = {lower=0, upper=0}
   for _, b in ipairs(self.boundary) do
      self.phiSheath[b]:clear(0.0)
      self.m0IonSheath[b]:clear(0.0)

      self.phiSheathPtr[b]   = self.phiSheath[b]:get(1)
      self.m0IonSheathPtr[b] = self.m0IonSheath[b]:get(1)

      self.boundaryIdxr[b] = self.phiSheath[b]:genIndexer()
   end

   if GKYL_HAVE_MPI then
      -- Need a communicator to broadcast the sheath potential and density along z.
      local commSet   = self.grid:commSet()
      local worldComm = commSet.comm
      local nodeComm  = commSet.nodeComm
      local nodeRank  = Mpi.Comm_rank(nodeComm)
      local zCommRank = 0
      if self.dim==3 then zCommRank = nodeRank % (self.grid:cuts(1)*self.grid:cuts(2)) end
      self.zComm = Mpi.Comm_split(worldComm, zCommRank, nodeRank)
      local zCommSize = Mpi.Comm_size(self.zComm)
      self.sheathRank = {lower=0, upper=zCommSize-1}
      self.numBoundaryDOFs = {lower=self.phiSheath["lower"]:localRange():volume()*self.phiSheath["lower"]:numComponents(),
                              upper=self.phiSheath["upper"]:localRange():volume()*self.phiSheath["upper"]:numComponents()}
      self.bcastSheathQuants = function(bInd)
         for d=1,self.dim do self.idxB[d] = 1 end
         self.phiSheath[bInd]:fill(self.boundaryIdxr[bInd](self.idxB), self.phiSheathPtr[bInd])
         self.m0IonSheath[bInd]:fill(self.boundaryIdxr[bInd](self.idxB), self.m0IonSheathPtr[bInd])
         Mpi.Bcast(self.phiSheathPtr[bInd]:data(), self.numBoundaryDOFs[bInd], Mpi.DOUBLE, self.sheathRank[bInd], self.zComm) 
         Mpi.Bcast(self.m0IonSheathPtr[bInd]:data(), self.numBoundaryDOFs[bInd], Mpi.DOUBLE, self.sheathRank[bInd], self.zComm) 
      end
   else
      self.bcastSheathQuants = function(bInd) end
   end
end

function ASheathPotential:_advance(tCurr, inFlds, outFlds)
   local GammaIon, m0Ion, jacobGeoInv = inFlds[1], inFlds[2], inFlds[3]
   local phi = outFlds[1]

   self.GammaIonPtr["lower"] = GammaIon["lower"]:get(1)
   self.GammaIonPtr["upper"] = GammaIon["upper"]:get(1)
   local m0IonPtr = m0Ion:get(1)
   local phiPtr   = phi:get(1)
   local jacobGeoInvPtr = jacobGeoInv:get(1)

   local innerIdxr = phi:genIndexer()

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
         self.idxB[self.sheathDir] = 1  -- Boundary grid has 1 cell in this direction.
         grid:getDx(self.dx)

         GammaIon[b]:fill(self.boundaryIdxr[b](self.idxB), self.GammaIonPtr[b])
         self.phiSheath[b]:fill(self.boundaryIdxr[b](self.idxB), self.phiSheathPtr[b])
         self.m0IonSheath[b]:fill(self.boundaryIdxr[b](self.idxB), self.m0IonSheathPtr[b])
         m0Ion:fill(innerIdxr(idx), m0IonPtr)
         jacobGeoInv:fill(innerIdxr(idx), jacobGeoInvPtr)
         
         self._phiSheathKer[b](self.dx[self.sheathDir], self.qElc, self.mElc, self.tempElc, self.halfSign[b], jacobGeoInvPtr:data(), self.GammaIonPtr[b]:data(), m0IonPtr:data(), self.m0IonSheathPtr[b]:data(), self.phiSheathPtr[b]:data())
      end

      -- Broadcast the sheath potential and density to other ranks along z.
      self.bcastSheathQuants(b)

      -- Loop over the half grid and compute phi.
      for idx in self.halfDomRangeDecomp[b]:rowMajorIter(tId) do
         idx:copyInto(self.idxB)
         self.idxB[self.sheathDir] = 1  -- Boundary grid has 1 cell in this direction.

         self.phiSheath[b]:fill(self.boundaryIdxr[b](self.idxB), self.phiSheathPtr[b])
         self.m0IonSheath[b]:fill(self.boundaryIdxr[b](self.idxB), self.m0IonSheathPtr[b])
         m0Ion:fill(innerIdxr(idx), m0IonPtr)
         phi:fill(innerIdxr(idx), phiPtr)
         jacobGeoInv:fill(innerIdxr(idx), jacobGeoInvPtr)

         self._phiKer(self.qElc, self.tempElc, jacobGeoInvPtr:data(), m0IonPtr:data(), self.m0IonSheathPtr[b]:data(), self.phiSheathPtr[b]:data(), phiPtr:data())
      end
   end
end

return ASheathPotential
