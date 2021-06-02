-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Gyrokinetic species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb, bcOpen)
-- and no additional setup.
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
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"

local GkBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkBasicBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   if tbl.kind=="function" or tbl.bcFunction then
      assert(type(tbl.bcFunction)=="function", "GkBasicBC: bcFunction must be a function.")
      self.bcKind   = "function"
      self.bcFuncIn = assert(tbl.bcFunction, "GkBasicBC: must specify the BC function in 'bcFunc' when using 'function' BC kind.")
      self.feedback = xsys.pickBool(tbl.feedback, false) 
   else
      self.bcKind = assert(tbl.kind, "GkBasicBC: must specify the type of BC in 'kind'.")

      self.phiWallFunc = tbl.phiWall
      if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "GkBasicBC: phiWall must be a function (t, xn).") end
   end
   self.evolve = xsys.pickBool(tbl.evolve, self.feedback or false) 

   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false
   if tbl.diagnostics then
     if #tbl.diagnostics>0 then self.saveFlux = true end
   end
end

function GkBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- Note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter
   for i = 1, self.basis:numBasis() do fOut[i] = 0.0 end
end
function GkBasicBC:bcOpen(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn, fOut)
end
function GkBasicBC:bcCopy(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   for i = 1, self.basis:numBasis() do fOut[i] = fIn[i] end
end
function GkBasicBC:bcReflect(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   self.basis:flipSign(dir, fIn, fOut)
   local vparDir = self.cdim+1
   self.basis:flipSign(vparDir, fOut, fOut)
end
function GkBasicBC:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f, fRefl)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.phiWallFld:fill(self.phiWallFldIdxr(idx), self.phiWallFldPtr)
   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_,
                                     self.phiPtr:data(), self.phiWallFldPtr:data(), f:data(), fRefl:data())
end
function GkBasicBC:bcSheath(dir, tm, idxIn, fIn, fOut)
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

function GkBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   -- Sheath BCs use phi and phiWall.
   self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   self.getPhi     = function(fieldIn, inIdx) return nil end 

   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   elseif self.bcKind == "open" then
      bcFunc   = function(...) return self:bcOpen(...) end
      skinType = "pointwise"
   elseif self.bcKind == "reflect" then
      assert(self.bcDir==self.cdim, "GkBasicBC: reflect BC can only be used along the last/parallel configuration space dimension.")
      bcFunc   = function(...) return self:bcReflect(...) end
      skinType = "flip"
   elseif self.bcKind == "function" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "sheath" then
      assert(self.bcDir==self.cdim, "GkBasicBC: sheath BC can only be used along the last/parallel configuration space dimension.")

      self.charge, self.mass = mySpecies.charge, mySpecies.mass

      self.fhatSheath = Lin.Vec(self.basis:numBasis())

      self.getPhi = function(fieldIn, inIdx) return fieldIn:rkStepperFields()[inIdx].phi end 
      -- Pre-create pointer and indexer for phi potential.
      local phi = field:rkStepperFields()[1].phi
      self.phiPtr, self.phiIdxr = phi:get(1), phi:genIndexer()

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
      if self.phiWallFunc then
         self.setPhiWall = Updater.EvalOnNodes {
            onGrid = field.grid,   evaluate = self.phiWallFunc,
            basis  = field.basis,  onGhosts = true,
         }
         self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
      end
      self.phiWallFld:sync(false)

      self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(self.basis:id(), self.cdim, 
                                                                             self.vdim, self.basis:polyOrder())

      bcFunc   = function(...) return self:bcSheath(...) end
      skinType = "flip"
   else
      assert(false, "GkBasicBC: BC kind not recognized.")
   end

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

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
      self:createBoundaryGrid()

      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid()

      local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
      -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
      self.allocCartField = function(self, grid, nComp, ghosts, metaData)
         local f = DataStruct.Field {
            onGrid        = grid,
            numComponents = nComp,
            ghost         = ghosts,
            metaData      = metaData,
         }
         f:clear(0.0)
         return f
      end
      local allocDistf = function()
         return self:allocCartField(self.boundaryGrid, self.basis:numBasis(),
                                    {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
      end
      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                    {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.basis:numBasis(),
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
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
         phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
         moment     = "GkM0", -- GkM0 = < f >
      }
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
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function GkBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkBasicBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with GkBasicBC.
function GkBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function GkBasicBC:getFlucF() return self.boundaryFluxRate end

function GkBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   self.setPhiWall:advance(tCurr, {}, self.phiWallFld) -- Compute wall potential if needed (i.e. sheath BC).
   self.phi = self.getPhi(field, inIdx)              -- If needed get the current plasma potential (for sheath BC).

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self.bcSolver:advance(tCurr, {fIn}, {fIn})
end

function GkBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local GkAbsorbBC = Proto(GkBasicBC)
function GkAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   GkAbsorbBC.super.fullInit(self, mySpecies)
end

local GkReflectBC = Proto(GkBasicBC)
function GkReflectBC:fullInit(mySpecies)
   self.tbl.kind  = "reflect"
   GkReflectBC.super.fullInit(self, mySpecies)
end

local GkCopyBC = Proto(GkBasicBC)
function GkCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   GkCopyBC.super.fullInit(self, mySpecies)
end

local GkOpenBC = Proto(GkBasicBC)
function GkOpenBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   GkOpenBC.super.fullInit(self, mySpecies)
end

local GkZeroFluxBC = Proto()
function GkZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end

local GkSheathBC = Proto(GkBasicBC)
function GkSheathBC:fullInit(mySpecies)
   self.tbl.kind  = "sheath"
   GkSheathBC.super.fullInit(self, mySpecies)
end
-- ................... End of GkBasicBC alias classes .................... --

return {GkBasic    = GkBasicBC,
        GkAbsorb   = GkAbsorbBC,
        GkCopy     = GkCopyBC,
        GkOpen     = GkOpenBC,
        GkReflect  = GkReflectBC,
        GkZeroFlux = GkZeroFluxBC,
        GkSheath   = GkSheathBC}
