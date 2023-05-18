-- Gkyl ------------------------------------------------------------------------
--
-- Boundary condition base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local Mpi        = require "Comm.Mpi"
local CartDecomp = require "Lib.CartDecomp"
local Grid       = require "Grid"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"

-- Empty shell source base class.
local BCsBase = Proto()

-- Functions that must be defined by subclasses.
function BCsBase:init(tbl) self.tbl = tbl end
function BCsBase:fullInit(speciesTbl) end
function BCsBase:setSpeciesName(nm) self.speciesName = nm end
function BCsBase:setName(nm) self.name = nm end
function BCsBase:setConfBasis(basis) self.confBasis = basis end
function BCsBase:setConfGrid(grid) self.confGrid = grid end
function BCsBase:setDir(dir) self.bcDir = dir end
function BCsBase:setEdge(edge) self.bcEdge = edge end
function BCsBase:setSaveFlux(newSaveFlux) self.saveFlux = newSaveFlux end
function BCsBase:bcCopy(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.basis:numBasis() do fOut[i] = fIn[i] end
end
function BCsBase:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. If setting them to exactly zero doesn't work it's likely that the
   -- sound speed & drift velocity diverge, so set them to something small.
   for i = 1, self.basis:numBasis() do fOut[i] = 0. end
end
function BCsBase:initCrossSpeciesCoupling(species) end
function BCsBase:createSolver(thisSpecies, extField) end
function BCsBase:createCouplingSolver(species, field, externalField) end
function BCsBase:createBoundaryGrid(ghostRange, ghostVec)
   -- Create a ghost boundary grid with only one cell in the direction of the BC.
   local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
   local reducedLowerRng, reducedUpperRng = {}, {}
   for d = 1, self.grid:ndim() do
      if d==self.bcDir then
         table.insert(reducedLower, self.bcEdge=="lower" and self.grid:lower(d)-ghostVec[d]*self.grid:dx(d)
                                                          or self.grid:upper(d))
         table.insert(reducedUpper, self.bcEdge=="lower" and self.grid:lower(d)
                                                          or self.grid:upper(d)+ghostVec[d]*self.grid:dx(d))
         table.insert(reducedNumCells, ghostVec[d])
         table.insert(reducedCuts, 1)
      else
         table.insert(reducedLower,    self.grid:lower(d))
         table.insert(reducedUpper,    self.grid:upper(d))
         table.insert(reducedNumCells, self.grid:numCells(d))
         table.insert(reducedCuts,     self.grid:cuts(d))
      end
      table.insert(reducedLowerRng, ghostRange:lower(d))
      table.insert(reducedUpperRng, ghostRange:upper(d))
   end
   local worldComm = self.grid:commSet().comm
   local worldRank = Mpi.Comm_rank(worldComm)
   local dirRank   = worldRank
   local cuts      = {}
   for d=1,3 do cuts[d] = self.grid:cuts(d) or 1 end
   local writeRank = Mpi.PROC_NULL
   if self.bcDir == 1 then
      dirRank = worldRank % (cuts[1]*cuts[2]) % cuts[1]
   elseif self.bcDir == 2 then
      dirRank = math.floor(worldRank/cuts[1]) % cuts[2]
   elseif self.bcDir == 3 then
      dirRank = math.floor(worldRank/cuts[1]/cuts[2])
   end
   self._splitComm = Mpi.Comm_split(worldComm, dirRank, worldRank)
   -- Set up which ranks to write from.
   if (self.bcEdge == "lower" and dirRank == 0) or
      (self.bcEdge == "upper" and dirRank == self.grid:cuts(self.bcDir)-1) then
      local worldGrp, splitGrp = Mpi.Comm_group(worldComm), Mpi.Comm_group(self._splitComm)
      local worldRanks = Lin.IntVec(1);  worldRanks[1] = worldRank
      writeRank = Mpi.Group_translate_ranks(worldGrp, worldRanks, splitGrp)[1]
   end
   self.writeRank = writeRank
   local reducedDecomp = CartDecomp.CartProd {
      comm      = self._splitComm,  cuts = reducedCuts,
      writeRank = writeRank,
   }
   self.boundaryGrid = Grid.RectCart {
      lower      = reducedLower,     cells         = reducedNumCells,
      upper      = reducedUpper,     decomposition = reducedDecomp,
      rangeLower = reducedLowerRng,  rangeUpper    = reducedUpperRng,
   }
end
function BCsBase:createConfBoundaryGrid(ghostRange, ghostVec)
   -- Create reduced boundary config-space grid with 1 cell in dimension of self._dir.
   if self.cdim == self.grid:ndim() then
      self.confBoundaryGrid = self.boundaryGrid
   else
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      local reducedLowerRng, reducedUpperRng = {}, {}
      for d = 1, self.cdim do
         if d==self.bcDir then
            table.insert(reducedLower, self.bcEdge=="lower" and self.grid:lower(d)-ghostVec[d]*self.grid:dx(d)
                                                             or self.grid:upper(d))
            table.insert(reducedUpper, self.bcEdge=="lower" and self.grid:lower(d)
                                                             or self.grid:upper(d)+ghostVec[d]*self.grid:dx(d))
            table.insert(reducedNumCells, ghostVec[d])
            table.insert(reducedCuts, 1)
         else
            table.insert(reducedLower,    self.grid:lower(d))
            table.insert(reducedUpper,    self.grid:upper(d))
            table.insert(reducedNumCells, self.grid:numCells(d))
            table.insert(reducedCuts,     self.grid:cuts(d))
         end
         table.insert(reducedLowerRng, ghostRange:lower(d))
         table.insert(reducedUpperRng, ghostRange:upper(d))
      end
      local reducedDecomp = CartDecomp.CartProd {
         comm      = self._splitComm,  cuts = reducedCuts,
         writeRank = self.writeRank,
      }
      self.confBoundaryGrid = Grid.RectCart {
         lower      = reducedLower,     cells         = reducedNumCells,
         upper      = reducedUpper,     decomposition = reducedDecomp,
         rangeLower = reducedLowerRng,  rangeUpper    = reducedUpperRng,
      }
   end
end
function BCsBase:createDiagnostics(mySpecies) end
function BCsBase:storeBoundaryFlux(tCurr, rkIdx, qOut) end
function BCsBase:calcCouplingMoments(tCurr, rkIdx, species) end
function BCsBase:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
function BCsBase:copyBoundaryFluxField(inIdx, outIdx) end
function BCsBase:combineBoundaryFluxField(outIdx, a, aIdx, ...) end
function BCsBase:computeBoundaryFluxRate(dtIn) end
function BCsBase:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx) end
function BCsBase:storeBoundaryFlux(tCurr, rkIdx, qOut) end
function BCsBase:copyBoundaryFluxField(inIdx, outIdx) end
function BCsBase:combineBoundaryFluxField(outIdx, a, aIdx, ...) end
function BCsBase:getDir() return self.bcDir end
function BCsBase:getEdge() return self.bcEdge end
function BCsBase:getKind() return self.bcKind end
function BCsBase:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self.bcEdge == "lower" then
      uv[self.bcDir] = global:lower(self.bcDir)-1   -- For ghost cells on "left".
   else
      lv[self.bcDir] = global:upper(self.bcDir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end
function BCsBase:evalOnBoundary(inFld)
   -- Evaluate inFld on the boundary grid and copy it to the boundary field buffer.
   self.boundaryField:copyRangeToRange(inFld, self.boundaryField:localRange(), self.myGlobalGhostRange)
   return self.boundaryField
end
function BCsBase:evalOnConfBoundary(inFld)
   -- For kinetic species this method evaluates inFld on the confBoundary grid.
   self.confBoundaryField:copyRangeToRange(inFld, self.confBoundaryField:localRange(), self.myGlobalConfGhostRange)
   return self.confBoundaryField
end

function BCsBase:createBoundaryTools(mySpecies,field,externalField) end

function BCsBase:createBoundaryToolsGK(mySpecies,field,externalField)
   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
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
   if self.saveflux or self.needsBoundaryTools then
      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
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
         self.boundaryFluxFields[i] = self.allocDistf()
      end
      self.boundaryFluxRate      = self.allocDistf()
      self.boundaryFluxFieldPrev = self.allocDistf()

      -- Part of global ghost range this rank owns.
      self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                      or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField = self:allocMoment()
      -- Range spanning ghost cells.
      self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                          or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
      local bmag = externalField.geo.bmag
      self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmag:getMetaData())
      self.bmag:copy(self:evalOnConfBoundary(bmag))
      local bmagInvSq = externalField.geo.bmagInvSq
      self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmagInvSq:getMetaData())
      self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq))
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
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
      self.phaseWeakMultiply = Updater.CartFieldBinOp {
         onGrid    = self.boundaryGrid,  operation = "Multiply",
         weakBasis = self.basis,         fieldBasis=self.confBasis,
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
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function BCsBase:createBoundaryToolsVlasov(mySpecies,field,externalField)
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


return BCsBase
