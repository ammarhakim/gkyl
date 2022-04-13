-- Gkyl ------------------------------------------------------------------------
--
-- Logical sheath boundary condition for a Gyrokinetic species.
-- See E. Shi, et al. PoP 22, 022504 (2015).
--
-- NOTE: MF 04/08/2022: currently only meant for 1x1v or 1x2v simulations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase      = require "App.BCs.BCsBase"
local DataStruct   = require "DataStruct"
local Updater      = require "Updater"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local Time         = require "Lib.Time"
local Range        = require "Lib.Range"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local CartDecomp   = require "Lib.CartDecomp"
local Grid         = require "Grid"
local DiagsApp     = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags      = require "App.Diagnostics.GkDiagnostics"
local xsys         = require "xsys"
local lume         = require "Lib.lume"
local root         = require "sci.root"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local DistFuncMomentCalcDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"

local GkLogicalSheathBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkLogicalSheathBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLogicalSheathBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind = "sheath"

   self.phiWallFunc = tbl.phiWall
   if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "GkLogicalSheathBC: phiWall must be a function (t, xn).") end

   self.evolve = xsys.pickBool(tbl.evolve, self.feedback or false) 

   self.saveFlux = true
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
      end
   end
end

function GkLogicalSheathBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkLogicalSheathBC:bcReflect(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   self.basis:flipSign(dir, fIn, fOut)
   local vparDir = self.cdim+1
   self.basis:flipSign(vparDir, fOut, fOut)
end
function GkLogicalSheathBC:cellAvTo0thDGcoeff(val, dgCoeffs)
   dgCoeffs[1] = val*math.sqrt(2)
end
function GkLogicalSheathBC:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f, fRefl)
   -- Set the sheath potential to phi_s = DeltaPhi + phi_wall = 0.5*m*vcut^2/|q|+phi_wall.
   self.phiWallFld:fill(self.phiWallFldIdxr(idx), self.phiWallFldPtr)
   local phiWallAtBcEdge = 0.
   for k = 1,self.confBasis:numBasis() do
      phiWallAtBcEdge = phiWallAtBcEdge
         + self.phiWallFldPtr[k]*self.confBasisAtBcEdgeSkin[k]
   end
   local phi_s = self.DeltaPhi + phiWallAtBcEdge
   self:cellAvTo0thDGcoeff(phi_s, self.phiPtr)

   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_,
                                     self.phiPtr:data(), self.phiWallFldPtr:data(), f:data(), fRefl:data())
end
function GkLogicalSheathBC:bcSheath(dir, tm, idxIn, fIn, fOut)
   -- Note that GK reflection only valid in z-vpar.
   -- This is checked when bc is created.

   -- Need to figure out if we are on lower or upper domain edge
   local edgeVal
   local globalRange = self.grid:globalRange()
   if idxIn[dir] == globalRange:lower(dir) then
      -- This means we are at lower domain edge,
      -- so we need to evaluate basis functions at z=-1.
      edgeVal = -1
   else
      -- This means we are at upper domain edge
      -- so we need to evaluate basis functions at z=1.
      edgeVal = 1
   end
   -- Get vpar limits of cell.
   local vpardir = self.cdim+1
   local gridIn  = self.grid
   gridIn:setIndex(idxIn)
   local vL = gridIn:cellLowerInDir(vpardir)
   local vR = gridIn:cellUpperInDir(vpardir)
   local vlowerSq, vupperSq
   -- This makes it so that we only need to deal with absolute values of vpar.
   if math.abs(vR)>=math.abs(vL) then
      vlowerSq = vL*vL
      vupperSq = vR*vR
   else
      vlowerSq = vR*vR
      vupperSq = vL*vL
   end
   local w  = gridIn:cellCenterInDir(vpardir)
   local dv = gridIn:dx(vpardir)
   -- calculate reflected distribution function fhat
   -- note: reflected distribution can be
   -- 1) fhat=0 (no reflection, i.e. absorb),
   -- 2) fhat=f (full reflection)
   -- 3) fhat=c*f (partial reflection)
   self:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, self.charge, self.mass, idxIn, fIn, self.fhatSheath)
   -- reflect fhat into skin cells
   self:bcReflect(dir, tm, nil, self.fhatSheath, fOut)
end

function GkLogicalSheathBC:initCrossSpeciesCoupling(species)
   for nm, s in lume.orderedIter(species) do
      if nm ~= self.speciesName then self.otherSpeciesName = nm end
   end

   -- Need to save a pointer to the other species' BC. Used later to fetch its boundary particle flux.
   local otherSpecies = species[self.otherSpeciesName]
   self.otherSpeciesBC = otherSpecies.nonPeriodicBCs[string.gsub(self.name,self.speciesName.."_","")]
end

function GkLogicalSheathBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local bcFunc, skinType
   assert(self.bcDir==self.cdim, "GkLogicalSheathBC: sheath BC can only be used along the last/parallel configuration space dimension.")

   self.charge, self.mass = mySpecies.charge, mySpecies.mass

   self.fhatSheath = Lin.Vec(self.basis:numBasis())

   -- Pre-create array holding phi in the skin cell.
   local numB = self.confBasis:numBasis()
   self.phiPtr = Lin.Vec(numB)
   for k = 1, numB do self.phiPtr[k] = 0. end

   -- Create field and function for calculating wall potential according to user-provided function.
   self.phiWallFld = DataStruct.Field {
      onGrid        = field.grid,
      numComponents = field.basis:numBasis(),
      ghost         = {1,1},
      metaData      = {polyOrder = field.basis:polyOrder(),
                       basisType = field.basis:id()},
      syncPeriodicDirs = false,
   }
   self.phiWallFld:clear(0.)
   self.phiWallFldPtr, self.phiWallFldIdxr = self.phiWallFld:get(1), self.phiWallFld:genIndexer()
   self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   if self.phiWallFunc then
      self.setPhiWall = Updater.EvalOnNodes {
         onGrid = field.grid,   evaluate = self.phiWallFunc,
         basis  = field.basis,  onGhosts = true,
      }
      self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
   end
   self.phiWallFld:sync(false)

   self.DeltaPhi = 0.

   self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(self.basis:id(), self.cdim, 
                                                                          self.vdim, self.basis:polyOrder())

   bcFunc   = function(...) return self:bcSheath(...) end
   skinType = "flip"

   local vdir = nil
   if self.bcDir==self.cdim then vdir = self.cdim+1 end

   self.bcSolver = Updater.Bc{
      onGrid   = self.grid,   edge               = self.bcEdge,  
      cdim     = self.cdim,   boundaryConditions = {bcFunc},   
      dir      = self.bcDir,  evaluate           = self.bcFuncIn,
      vdir     = vdir,        evolveFn           = self.evolve,
      skinLoop = skinType,    feedback           = self.feedback,
      basis    = self.basis,  confBasis          = self.confBasis,
      advanceArgs = {{mySpecies:rkStepperFields()[1]}, {mySpecies:rkStepperFields()[1]}},
   }

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid(mySpecies)

   -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
   self:createConfBoundaryGrid(mySpecies)

   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
   self.allocMoment = function(self)
      return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                 {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
   end

   self.fluxM0 = {self = self:allocMoment()}  -- This species M0 moment of the boundary flux, used to find vcut.
   self.fluxM0Ptr, self.fluxM0Idxr  = {self = self.fluxM0["self"]:get(1)}, {self = self.fluxM0["self"]:genIndexer()}
   self.fluxM0atBcEdge = {self = 0., other = 0.}

   self.fluxM0partial = Lin.Vec(self.confBasis:numBasis())  -- Partial M0 moment of the boundary flux.
   self.fluxM0partialPrev = Lin.Vec(self.confBasis:numBasis())
   self.xcP, self.dxP = Lin.Vec(self.ndim), Lin.Vec(self.ndim)

   -- Quantities used in evaluating the particle flux at the boundary.
   self.confBasisAtBcEdgeGhost, self.confBasisAtBcEdgeSkin = Lin.Vec(numB), Lin.Vec(numB)
   self.confBasis:evalBasis({self.bcEdge=="lower" and 1. or -1.}, self.confBasisAtBcEdgeGhost)
   self.confBasis:evalBasis({self.bcEdge=="lower" and -1. or 1.}, self.confBasisAtBcEdgeSkin)

   -- Set parameters needed by cut-off velocity (vcut) calculation.
   local phaseRange = distf:localRange()
   -- Lower/upper limits of velocity space loops when searching for vcut.
   -- Here we assume the vpar grid has even number of cells and is symmetric about vpar=0.
   self.vLoIdx, self.vUpIdx, self.vStep = {}, {}, {}
   self.vUpIdxPartial = {}
   if self.bcEdge == "lower" then
      self.vLoIdx[1], self.vUpIdx[1] = phaseRange:lower(self.cdim+1), phaseRange:shape(self.cdim+1)/2
      self.vStep[1] = 1
      self.set_vUpIdxPartial = function(cutIdx) self.vUpIdxPartial[1] = cutIdx-1 end
   elseif self.bcEdge == "upper" then
      self.vLoIdx[1], self.vUpIdx[1] = phaseRange:upper(self.cdim+1), phaseRange:shape(self.cdim+1)/2+1
      self.vStep[1] = -1
      self.set_vUpIdxPartial = function(cutIdx) self.vUpIdxPartial[1] = cutIdx+1 end
   end
   if self.vdim == 1 then
      self.vLoIdx[2], self.vUpIdx[2] = 1, 1
      self.vStep[2] = 1
   elseif self.vdim == 2 then
      self.vLoIdx[2], self.vUpIdx[2] = phaseRange:lower(self.cdim+2), phaseRange:upper(self.cdim+2)
      self.vStep[2] = 1
   end

   self._m0Ker        = DistFuncMomentCalcDecl.selectGkMomCalc("GkM0", self.confBasis:id(), self.cdim,
                                                                self.vdim, self.basis:polyOrder())
   self._m0partialKer = DistFuncMomentCalcDecl.selectGkMomPartialCalc("GkM0", self.confBasis:id(), self.cdim,
                                                                       self.vdim, self.basis:polyOrder())

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

      local allocDistf = function()
         return self:allocCartField(self.boundaryGrid, self.basis:numBasis(),
                                    {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(),
                                    {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.boundaryPtr        = {}
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
         self.boundaryPtr[i]        = self.boundaryFluxFields[i]:get(1)
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()
      self.boundaryIdxr          = self.boundaryFluxFields[1]:genIndexer()

      self.idxOut = Lin.IntVec(self.grid:ndim())

      -- Create the range needed to loop over ghosts.
      local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
      self.ghostRange = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRange, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField    = self:allocMoment()
      self.confBoundaryFieldPtr = self.confBoundaryField:get(1)
      self.confBoundaryIdxr     = self.confBoundaryField:genIndexer()
      local confGlobal        = numDensity:globalRange()
      local confGlobalExt     = numDensity:globalExtRange()
      local confLocalExtRange = numDensity:localExtRange()
      self.confGhostRange = confLocalExtRange:intersect(self:getGhostRange(confGlobal, confGlobalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self.confGhostRangeDecomp = LinearDecomp.LinearDecompRange {range=self.confGhostRange, numSplit=self.grid:numSharedProcs()}

      -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
      local bmag = externalField.geo.bmag 
      self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                      {bmag:lowerGhost(),bmag:upperGhost()}, bmag:getMetaData())
      self.bmag:copy(self:evalOnConfBoundary(bmag))
      self.bmagPtr = self.bmag:get(1)
      local bmagInvSq = externalField.geo.bmagInvSq
      self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                          {bmagInvSq:lowerGhost(),bmagInvSq:upperGhost()}, bmagInvSq:getMetaData())
      self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq))
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                             {jacobGeo:lowerGhost(),jacobGeo:upperGhost()}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                                {jacobGeoInv:lowerGhost(),jacobGeoInv:upperGhost()}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

      -- Declare methods/functions needed for handling saved fluxes and needed by diagnostics.
      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRangeDecomp:rowMajorIter(self.tId) do
            idx:copyInto(self.idxOut)
            qOut:fill(self.distfInIdxr(idx), ptrOut)

            -- Before operating on ghosts, store ghost values for later flux diagnostics
            self.idxOut[self.bcDir] = 1
            self.boundaryFluxFields[rkIdx]:fill(self.boundaryIdxr(self.idxOut), self.boundaryPtr[rkIdx])
            for c = 1, qOut:numComponents() do self.boundaryPtr[rkIdx][c] = ptrOut[c] end
         end
      end
      self.copyBoundaryFluxFieldFunc = function(inIdx, outIdx)
         self.boundaryFluxFields[outIdx]:copy(self.boundaryFluxFields[inIdx])
      end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...)
         local args  = {...} -- Package up rest of args as table.
         local nFlds = #args/2
         self.boundaryFluxFields[outIdx]:combine(a, self.boundaryFluxFields[aIdx])
         for i = 1, nFlds do -- Accumulate rest of the fields.
            self.boundaryFluxFields[outIdx]:accumulate(args[2*i-1], self.boundaryFluxFields[args[2*i]])
         end
      end

      -- Number density calculator. Needed regardless of diagnostics.
      local mass = mySpecies.mass
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
         phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
         moment     = "GkM0", -- GkM0 = < f >
      }

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         self.calcBoundaryFluxRateFunc = function(dtIn)
            -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
            self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                          -1.0/dtIn, self.boundaryFluxFieldPrev)
            self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
         end
         -- Set up weak multiplication and division operators (for diagnostics).
         self.confWeakMultiply = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Multiply",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Divide",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,         quantity      = "V",
            }
         }
         -- Moment calculators (for diagnostics).
         local mass = mySpecies.mass
         self.momDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM1", -- GkM1 = < v_parallel f >
         }
         self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM2", -- GkM2 = < (v_parallel^2 + 2*mu*B/m) f >
         }
         self.M2parCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM2par", -- GkM2par = < v_parallel^2 f >
         }
         self.M3parCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM3par", -- GkM3par = < v_parallel^3 f >
         }
         if self.vdim > 1 then
            self.M2perpCalc = Updater.DistFuncMomentCalc {
               onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
               phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
               moment     = "GkM2perp", -- GkM2 = < (mu*B/m) f >
            }
            self.M3perpCalc = Updater.DistFuncMomentCalc {
               onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
               phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
               moment     = "GkM3perp", -- GkM3perp = < vpar*(mu*B/m) f >
            }
         end
         self.divideByJacobGeo = self.jacobGeoInv
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
         self.multiplyByJacobGeo = self.jacobGeo
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function GkLogicalSheathBC:createCouplingSolver(species, field, externalField)
   -- Fetch and organize information about other species that this BCs depends on.
   local otherSpecies = species[self.otherSpeciesName]

   -- Need to save a pointer to the other species' BC. Used later to fetch its boundary particle flux.
   self.otherSpeciesBC = otherSpecies.nonPeriodicBCs[string.gsub(self.name,self.speciesName.."_","")]
   self.fluxM0["other"] = self.otherSpeciesBC.fluxM0["self"]
   self.fluxM0Ptr["other"], self.fluxM0Idxr["other"]  = self.fluxM0["other"]:get(1), self.fluxM0["other"]:genIndexer()
end

function GkLogicalSheathBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkLogicalSheathBC:evalBoundaryConfFieldAtBcEdge(fPtrIn)
   -- Given a conf-space field 'fIn' defined confBoundaryGrid, with the pointer 'fPtrIn' to
   -- the field in a particular cell, evaluate such field at the lower/upper boundary
   -- of the cell depending on the value of bcEdge.
   -- Note: recall that the boundary flux needs to be multiplied by dx/2 in the direction of the BC.
   self.boundaryGrid:getDx(self.dxP)
   local fInAtBcEdge = 0.
   for k = 1,self.confBasis:numBasis() do
      fInAtBcEdge = fInAtBcEdge
         + (0.5*self.dxP[self.cdim])*fPtrIn[k]*self.confBasisAtBcEdgeGhost[k]
   end
   return fInAtBcEdge
end

function GkLogicalSheathBC:m0Integral_at_vparCell(jCurr, fFlux, fFluxPtr)
   for k = self.vLoIdx[2], self.vUpIdx[2] do

      local idx = {1,jCurr,k}
      self.boundaryGrid:setIndex(idx)
      self.boundaryGrid:getDx(self.dxP)
      self.boundaryGrid:cellCenter(self.xcP)

      self.bmag:fill(self.confBoundaryIdxr(idx), self.bmagPtr)
      fFlux:fill(self.boundaryIdxr(idx), fFluxPtr)

      self._m0Ker(self.xcP:data(), self.dxP:data(), self.mass, self.bmagPtr:data(), fFluxPtr:data(), self.fluxM0partial:data())
   end
end

-- Stopping criteria for root-finding.
local function stop(tol)
   return function(x, y, xl, xu, yl, yu)
      if math.abs(y) < tol then return true else return false end
   end
end

function GkLogicalSheathBC:calcCouplingMoments(tCurr, rkIdx, species) end

function GkLogicalSheathBC:calcDeltaPhi(rkIdx)
   -- Compute the change in the potential across the sheath:
   --   DeltaPhi = phi_s - phi_wall = 0.5*m*vcut^2/q - phi_wall

   -- Relative tolerance (relative to the particle flux of the absorbed species)
   -- used when trying to find the cutoff velocity.
   local relTol = 1.e-11

   -- Establish the partially reflected species:
   for _, sp in ipairs({"self","other"}) do
      self.fluxM0[sp]:fill(self.fluxM0Idxr[sp]({1}), self.fluxM0Ptr[sp])
      self.fluxM0atBcEdge[sp] = self:evalBoundaryConfFieldAtBcEdge(self.fluxM0Ptr[sp])
   end
--   print(string.format("%s %s fluxes: %s = %g | %s = %g | DeltaPhi=%g",self.name, self.speciesName, self.speciesName, self.fluxM0atBcEdge["self"],self.otherSpeciesName,self.fluxM0atBcEdge["other"],self.DeltaPhi))

   if self.fluxM0atBcEdge["self"] > self.fluxM0atBcEdge["other"] then
      -- Partially reflect this species.

      -- Find the cell containing the cutoff velocity.
      local vcut, vcutIdx = nil, nil
      local fFlux    = self:getBoundaryFluxFields()[rkIdx]
      local fFluxPtr = fFlux:get(1)
      for k = 1,self.confBasis:numBasis() do self.fluxM0partial[k] = 0. end
--      print("at first search ",self.speciesName,self.vLoIdx[1], self.vUpIdx[1], self.vStep[1])
      for j = self.vLoIdx[1], self.vUpIdx[1], self.vStep[1] do

         self:m0Integral_at_vparCell(j, fFlux, fFluxPtr)

         -- Evaluate particle flux at the boundary:
         local fluxM0partialAtBcEdge = self:evalBoundaryConfFieldAtBcEdge(self.fluxM0partial)

         if math.abs(fluxM0partialAtBcEdge-self.fluxM0atBcEdge["other"]) < relTol*self.fluxM0atBcEdge["other"] then
            -- The cutoff velocity is the upper boundary along vpar of this cell.
            vcut = self.xcP[2]+0.5*self.dxP[2]
         elseif fluxM0partialAtBcEdge > self.fluxM0atBcEdge["other"] then
            -- We just passed vcut, so this is the cell containing vcut.
            vcutIdx = j
            break
         elseif j == self.vUpIdx[1] then
            -- The cutoff velocity is above the maximum vpar. Set it to something slightly above.
--            print("end ",self.name, j, " partialFlux=",fluxM0partialAtBcEdge)
--            print("end ",self.name, j, " partialFlux=",self.fluxM0partial[1], self.fluxM0partial[2])
            vcut = self.bcEdge=="lower" and 1.5*self.boundaryGrid:lower(2) or 1.5*self.boundaryGrid:upper(2)
         end
      end
   
      if vcutIdx then
         -- vcut is inside the vcutIdx cell (most common case).
         -- Find it with a root finding method.

         -- Integrate up to vcutIdx-1 if bcEdge=lower or vcutIdx+1 if bcEdge=upper:
         self.set_vUpIdxPartial(vcutIdx)
         for k = 1,self.confBasis:numBasis() do self.fluxM0partial[k] = 0. end
         for j = self.vLoIdx[1], self.vUpIdxPartial[1], self.vStep[1] do
            self:m0Integral_at_vparCell(j, fFlux, fFluxPtr)
         end

         local function rootEq(vcutIn)
            for k = 1,self.confBasis:numBasis() do self.fluxM0partialPrev[k] = self.fluxM0partial[k] end
            -- The contribution to self.fluxM0partial from this cell in vpar:
            for k = self.vLoIdx[2], self.vUpIdx[2] do
               local idx = {1,vcutIdx,k}
               self.boundaryGrid:setIndex(idx)
               self.boundaryGrid:getDx(self.dxP)
               self.boundaryGrid:cellCenter(self.xcP)

               self.bmag:fill(self.confBoundaryIdxr(idx), self.bmagPtr)
               fFlux:fill(self.boundaryIdxr(idx), fFluxPtr)

               self._m0partialKer[self.bcEdge](self.xcP:data(), self.dxP:data(), self.mass, self.bmagPtr:data(),
                                               vcutIn, fFluxPtr:data(), self.fluxM0partialPrev:data())
            end

            -- Evaluate the particle flux at the boundary
            local fluxM0partialAtBcEdge = self:evalBoundaryConfFieldAtBcEdge(self.fluxM0partialPrev)

            return fluxM0partialAtBcEdge-self.fluxM0atBcEdge["other"]
         end

         local idx_vpar = {1,vcutIdx,1}
         self.boundaryGrid:setIndex(idx_vpar)
         self.boundaryGrid:getDx(self.dxP)
         self.boundaryGrid:cellCenter(self.xcP)
         local vpar_c, DvparD2 = self.xcP[self.cdim+1], self.dxP[self.cdim+1]/2.

         vcut = root.ridders(rootEq, vpar_c-DvparD2, vpar_c+DvparD2, stop(relTol*self.fluxM0atBcEdge["other"]))
      end
      
      -- Having found vcut comput DeltaPhi = phiSheath-phiWall = m*vcut^2/2
--      print(string.format("%s vcutIdx = %d | vcut = %g",self.name,vcutIdx==nil and -99 or vcutIdx,vcut))
      self.DeltaPhi = 0.5*self.mass*(vcut^2)/math.abs(self.charge)

      -- Set DeltaPhi in other other species, which will use it to apply
      -- sheath BCs but due to its charge it will not be reflected at all.
      self.otherSpeciesBC.DeltaPhi = self.DeltaPhi
   end
end

function GkLogicalSheathBC:advanceCrossSpeciesCoupling(tCurr, species, outIdx)
   -- Compute the 0th moment of the boundary flux.
   self.numDensityCalc:advance(tCurr, {self:getBoundaryFluxFields()[outIdx]}, {self.fluxM0["self"]})
   self.otherSpeciesBC.numDensityCalc:advance(tCurr, {self.otherSpeciesBC:getBoundaryFluxFields()[outIdx]}, {self.fluxM0["other"]})

   self:calcDeltaPhi(outIdx)  -- Compute change in phi across the sheath.
end

function GkLogicalSheathBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkLogicalSheathBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkLogicalSheathBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkLogicalSheathBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the GkDiagnostics with GkLogicalSheathBC.
function GkLogicalSheathBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function GkLogicalSheathBC:getFlucF() return self.boundaryFluxRate end

function GkLogicalSheathBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   self.setPhiWall:advance(tCurr, {}, {self.phiWallFld}) -- Compute wall potential if needed (i.e. sheath BC).

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self.bcSolver:advance(tCurr, {fIn}, {fIn})
end

function GkLogicalSheathBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return GkLogicalSheathBC
