local Proto      = require "Lib.Proto"
local BCsBase    = require "App.BCs.BCsBase"
local Updater    = require "Updater"
local DataStruct = require "DataStruct"

local BCsBaseVlasov = Proto(BCsBase)

function BCsBaseVlasov:createBoundaryToolsVlasov(mySpecies,field,externalField)
   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   local globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   self:createBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   self.allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, distf:getMetaData())
   end
   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux or self.needsBoundaryTools then
      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntThreeMoments = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.vdim+2, {0,0}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntMoment = function(self, comp)
         local metaData = {charge = self.charge,  mass = self.mass,}
         local ncomp = comp or 1
         local f = DataStruct.DynVector{numComponents = ncomp,     writeRank = self.confBoundaryGrid:commSet().writeRank,
                                        metaData      = metaData,  comm      = self.confBoundaryGrid:commSet().comm,}
         return f
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()

      -- Part of global ghost range this rank owns.
      self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                      or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField = self:allocMoment()
      -- Range spanning ghost cells.
      self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                          or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      local jacobGeo = externalField.geo and externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo and externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

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
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
         phaseBasis = self.basis,         moment     = "M0",
      }

      -- Integrated number density calculator. Needed regardless of diagnostics (for steady state sources).
      self.integNumDensityCalc = Updater.DistFuncMomentDG {
         onGrid     = self.boundaryGrid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,          moment     = "M0",
         isIntegrated = true,
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
         self.confWeakDotProduct = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "DotProduct",
            onGhosts  = true,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,         quantity      = "V",
            },
            vector = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = self.vdim,
               basis  = self.confBasis,         quantity      = "V",
            },
         }
         -- Moment calculators (for diagnostics).
         self.momDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M1i",
         }
         self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M2",
         }
         self.M2ijCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M2ij",
         }
         self.M3iCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M3i",
         }
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

return BCsBaseVlasov

