-- Gkyl ------------------------------------------------------------------------
--
-- Boundary condition base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local CartDecomp = require "Lib.CartDecomp"
local Grid       = require "Grid"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"
local Mpi        = require "Comm.Mpi"

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
         table.insert(reducedLowerRng, ghostRange:lower(d))
         table.insert(reducedUpperRng, ghostRange:upper(d))
      else
         table.insert(reducedLower,    self.grid:lower(d))
         table.insert(reducedUpper,    self.grid:upper(d))
         table.insert(reducedNumCells, self.grid:numCells(d))
         table.insert(reducedCuts,     self.grid:cuts(d))
         table.insert(reducedLowerRng, self.grid:globalRange():lower(d))
         table.insert(reducedUpperRng, self.grid:globalRange():upper(d))
      end
   end
   local parentComm = self.grid:commSet().host
   local parentRank = Mpi.Comm_rank(parentComm)
   local dirRank   = parentRank
   local cuts      = {}
   for d=1,3 do cuts[d] = self.grid:cuts(d) or 1 end
   local writeRank = Mpi.PROC_NULL
   if self.bcDir == 1 then
      dirRank = parentRank % (cuts[1]*cuts[2]) % cuts[1]
   elseif self.bcDir == 2 then
      dirRank = math.floor(parentRank/cuts[1]) % cuts[2]
   elseif self.bcDir == 3 then
      dirRank = math.floor(parentRank/cuts[1]/cuts[2])
   end

   local msn = self.grid:getMessenger()
   self._splitComm = Mpi.Comm_split(parentComm, dirRank, parentRank)
   if GKYL_HAVE_CUDA then
      self._splitComm_dev, self.ncclIdSplit = msn:newCommDevice(self._splitComm)
   end

   -- Set up which ranks to write from.
   if (self.bcEdge == "lower" and dirRank == 0) or
      (self.bcEdge == "upper" and dirRank == self.grid:cuts(self.bcDir)-1) then
      local parentGrp, splitGrp = Mpi.Comm_group(parentComm), Mpi.Comm_group(self._splitComm)
      local parentRanks = Lin.IntVec(1);  parentRanks[1] = parentRank
      writeRank = Mpi.Group_translate_ranks(parentGrp, parentRanks, splitGrp)[1]
   end
   self.writeRank = writeRank
   local reducedDecomp = CartDecomp.CartProd {
      cuts  = reducedCuts,  writeRank = writeRank,
      comms = {host=self._splitComm, device=self._splitComm_dev,},
   }
   self.boundaryGrid = Grid.RectCart {
      lower      = reducedLower,     cells         = reducedNumCells,
      upper      = reducedUpper,     decomposition = reducedDecomp,
      rangeLower = reducedLowerRng,  rangeUpper    = reducedUpperRng, 
      messenger  = msn,
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
            table.insert(reducedLowerRng, ghostRange:lower(d))
            table.insert(reducedUpperRng, ghostRange:upper(d))
         else
            table.insert(reducedLower,    self.grid:lower(d))
            table.insert(reducedUpper,    self.grid:upper(d))
            table.insert(reducedNumCells, self.grid:numCells(d))
            table.insert(reducedCuts,     self.grid:cuts(d))
            table.insert(reducedLowerRng, self.grid:globalRange():lower(d))
            table.insert(reducedUpperRng, self.grid:globalRange():upper(d))
         end
      end
      local reducedDecomp = CartDecomp.CartProd {
         cuts  = reducedCuts,  writeRank = self.writeRank,
         comms = {host=self._splitComm, device=self._splitComm_dev,},
      }
      self.confBoundaryGrid = Grid.RectCart {
         lower      = reducedLower,     cells         = reducedNumCells,
         upper      = reducedUpper,     decomposition = reducedDecomp,
         rangeLower = reducedLowerRng,  rangeUpper    = reducedUpperRng, 
         messenger  = self.grid:getMessenger(),
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
function BCsBase:evalOnBoundary(inFld, outFld)
   -- Evaluate inFld on the boundary grid and copy it to the boundary field buffer.
   outFld:copyRangeToRange(inFld, outFld:localRange(), self.myGlobalGhostRange)
   return outFld
end
function BCsBase:evalOnConfBoundary(inFld, outFld)
   -- For kinetic species this method evaluates inFld on the confBoundary grid.
   outFld:copyRangeToRange(inFld, outFld:localRange(), self.myGlobalConfGhostRange)
   return outFld
end

return BCsBase
