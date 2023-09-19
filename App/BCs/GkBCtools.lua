-- Gkyl ------------------------------------------------------------------------
--
-- BC related tools for gyrokinetics. Mostly to avoid code duplication
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local Updater    = require "Updater"

local _M = {}

function _M.createBoundaryTools(mySpecies, field, externalField, bcApp)
   local self = bcApp

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()

   local globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

   -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
   self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData, syncPeriodic)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
         syncPeriodicDirs = syncPeriodic,
      }
      f:clear(0.0)
      return f
   end
   self.allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0},
                                 {polyOrder = self.basis:polyOrder(), basisType = self.basis:id()})
   end
   self.allocMoment = function(self)
      return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
   end
   self.allocVectorMoment = function(self, dim)
      return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
   end
   self.allocIntMoment = function(self, comp)
      local metaData = {charge = self.charge,  mass = self.mass,}
      local ncomp = comp or 1
      local gridWriteRank = self.confBoundaryGrid:commSet().writeRank
      local f = DataStruct.DynVector{numComponents = ncomp,     writeRank = gridWriteRank<0 and gridWriteRank or 0,
                                     metaData      = metaData,  comm      = self.confBoundaryGrid:commSet().comm,}
      return f
   end

   -- Part of global ghost range this rank owns.
   self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                   or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

   -- Range spanning ghost cells.
   self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                       or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]
end

function _M.createFluxTools(mySpecies, field, externalField, bcApp)
   local self = bcApp

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()

   -- Allocate fields needed.
   self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
   self.distfInIdxr        = distf:genIndexer()
   for i = 1, #mySpecies:rkStepperFields() do
      self.boundaryFluxFields[i] = self.allocDistf()
   end
   self.boundaryFluxRate      = self.allocDistf()
   self.boundaryFluxFieldPrev = self.allocDistf()

   -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
   self.confBoundaryField = self:allocMoment()

   -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
   local bmag = externalField.geo.bmag
   self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmag:getMetaData())
   self.bmag:copy(self:evalOnConfBoundary(bmag, self.confBoundaryField))
   local bmagInvSq = externalField.geo.bmagInvSq
   self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmagInvSq:getMetaData())
   self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq, self.confBoundaryField))
   local jacobGeo = externalField.geo.jacobGeo
   if jacobGeo then
      self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
      self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo, self.confBoundaryField))
   end
   local jacobGeoInv = externalField.geo.jacobGeoInv
   if jacobGeoInv then
      self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
      self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv, self.confBoundaryField))
   end

   -- Declare methods/functions needed for handling saved fluxes and needed by diagnostics.
   self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
      self.boundaryFluxFields[rkIdx]:copyRangeToRange(qOut, self.boundaryFluxFields[rkIdx]:localRange(), self.myGlobalGhostRange)
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
end

function _M.createDiagnosticTools(mySpecies, field, externalField, bcApp)
   local self = bcApp

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

return _M
