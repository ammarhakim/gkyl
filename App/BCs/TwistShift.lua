-- Gkyl ------------------------------------------------------------------------
--
-- Apply twist shift BCs for flux-tube simulations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"
local Grid       = require "Grid"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"
local lume       = require "Lib.lume"

local TwistShiftBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function TwistShiftBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function TwistShiftBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.yShiftFuncIn = assert(tbl.shiftFunction, "TwistShiftBC: must provide the function that computes the y-shift in 'shiftFunction'.")
   
   self.yShiftPolyOrder = tbl.shiftPolyOrder

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function TwistShiftBC:setName(nm) self.name = self.speciesName.."_"..nm end

function TwistShiftBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   if self.yShiftPolyOrder == nil then self.yShiftPolyOrder = self.basis:polyOrder() end

   if self.bcEdge=="lower" then
      self.yShiftFunc = function(t,xn) return -self.yShiftFuncIn(t,xn) end
   else
      self.yShiftFunc = self.yShiftFuncIn
   end

   self.bcSolver = Updater.TwistShiftBC {
      onGrid    = self.grid,       yShiftFunc      = self.yShiftFunc,
      basis     = self.basis,      yShiftPolyOrder = self.yShiftPolyOrder,
      confBasis = self.confBasis,  edge            = self.bcEdge,
   }

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid()

   local distf = mySpecies["getDistF"] and mySpecies:getDistF() or mySpecies:getMoments()
   -- Define methods to allocate fields defined on boundary grid (used by e.g. diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(),
                                 {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
   end
   self.boundaryField    = allocDistf()
   self.boundaryFieldPtr = self.boundaryField:get(1)

   -- Create the range needed to loop over ghosts.
   local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
   self.ghostRange = localExtRange:intersect(self:getGhostRange(global, globalExt))

   self.boundaryIdxr = self.boundaryField:genIndexer()

   self.idxOut = Lin.IntVec(self.grid:ndim())

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid()

      local numDensity = mySpecies:getNumDensity()
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

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField    = self:allocMoment()
      self.confBoundaryFieldPtr = self.confBoundaryField:get(1)
      self.confBoundaryIdxr     = self.confBoundaryField:genIndexer()
      local confGlobal        = numDensity:globalRange()
      local confGlobalExt     = numDensity:globalExtRange()
      local confLocalExtRange = numDensity:localExtRange()
      self.confGhostRange = confLocalExtRange:intersect(self:getGhostRange(confGlobal, confGlobalExt)) -- Range spanning ghost cells.

      -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
      local bmag = externalField.geo.bmag 
      self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                      {bmag:lowerGhost(),bmag:upperGhost()}, bmag:getMetaData())
      self.bmag:copy(self:evalOnConfBoundary(bmag, self.confBoundaryField))
      local bmagInvSq = externalField.geo.bmagInvSq
      self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                          {bmagInvSq:lowerGhost(),bmagInvSq:upperGhost()}, bmagInvSq:getMetaData())
      self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq, self.confBoundaryField))
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                             {jacobGeo:lowerGhost(),jacobGeo:upperGhost()}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo, self.confBoundaryField))
      end
      local jacobGeoInv = externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                                {jacobGeoInv:lowerGhost(),jacobGeoInv:upperGhost()}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv, self.confBoundaryField))
      end

      -- Declare methods/functions needed for handling saved fluxes and needed by diagnostics.
      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRange:rowMajorIter() do
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

      -- Number density calculator. Needed regardless of diagnostics (for recycling BCs).
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
            weakBasis = self.confBasis,  operation = "Multiply",
            onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Divide",
            onRange   = self.confBoundaryField:localRange(),  onGhosts = false,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegrate {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,
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

function TwistShiftBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function TwistShiftBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function TwistShiftBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function TwistShiftBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function TwistShiftBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with TwistShiftBC.
function TwistShiftBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                self.boundaryFluxRate, self.boundaryFluxRate} end
function TwistShiftBC:getFlucF() return self.boundaryFluxRate end

function TwistShiftBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   -- Apply periodic BCs in z direction, but only on the first edge (lower or upper, whichever is first)
   -- First check if sync has been called on this step by any BC.
   local doSync = true
   for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do 
      if bc.synced then doSync = false end
   end
   if doSync then
      -- Fetch skin cells from opposite z-edge (donor) and copy to ghost cells (target).
      local fldGrid = fIn:grid()
      local periodicDirs = fldGrid:getPeriodicDirs()
      local modifiedPeriodicDirs = {3}
      fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
      fIn:sync()
      fldGrid:setPeriodicDirs(periodicDirs)
      self.synced = true
   else
      -- Don't sync but reset synced flags for next step.
      for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do bc.synced = false end
   end

   -- Copy field ghost cells to boundary field (buffer).
   self:evalOnBoundary(fIn, self.boundaryField)

   self.bcSolver:advance(tCurr, {self.boundaryField}, {fIn})
end

function TwistShiftBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return TwistShiftBC
