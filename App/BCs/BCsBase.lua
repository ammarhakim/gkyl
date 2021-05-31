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
function BCsBase:createSolver(thisSpecies, extField) end
function BCsBase:createBoundaryGrid()
   -- Create a ghost boundary grid with only one cell in the direction of the BC.
   if self.grid:isShared() then
      assert(false, "VlasovBasicBC: shared memory implementation of boundary flux diagnostics not ready.")
   else
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      for d = 1, self.grid:ndim() do
         if d==self.bcDir then
            table.insert(reducedLower, -self.grid:dx(d)/2.)
            table.insert(reducedUpper,  self.grid:dx(d)/2.)
            table.insert(reducedNumCells, 1)
            table.insert(reducedCuts, 1)
         else
            table.insert(reducedLower,    self.grid:lower(d))
            table.insert(reducedUpper,    self.grid:upper(d))
            table.insert(reducedNumCells, self.grid:numCells(d))
            table.insert(reducedCuts,     self.grid:cuts(d))
         end
      end
      local commSet  = self.grid:commSet()
      local worldComm, nodeComm = commSet.comm, commSet.nodeComm
      local nodeRank = Mpi.Comm_rank(nodeComm)
      local dirRank  = nodeRank
      local cuts     = {}
      for d=1,3 do cuts[d] = self.grid:cuts(d) or 1 end
      local writeRank = -1
      if self.bcDir == 1 then
         dirRank = nodeRank % (cuts[1]*cuts[2]) % cuts[1]
      elseif self.bcDir == 2 then
         dirRank = math.floor(nodeRank/cuts[1]) % cuts[2]
      elseif self.bcDir == 3 then
         dirRank = math.floor(nodeRank/cuts[1]/cuts[2])
      end
      self._splitComm = Mpi.Comm_split(worldComm, dirRank, nodeRank)
      -- Set up which ranks to write from.
      if self.bcEdge == "lower" and dirRank == 0 then
         writeRank = nodeRank
      elseif self.bcEdge == "upper" and dirRank == self.grid:cuts(self.bcDir)-1 then
         writeRank = nodeRank
      end
      self.writeRank = writeRank
      local reducedDecomp = CartDecomp.CartProd {
         comm      = self._splitComm,  cuts      = reducedCuts,
         writeRank = writeRank,        useShared = self.grid:isShared(),
      }
      self.boundaryGrid = Grid.RectCart {
         lower = reducedLower,  cells         = reducedNumCells,
         upper = reducedUpper,  decomposition = reducedDecomp,
      }
   end
end
function BCsBase:createConfBoundaryGrid()
-- Create reduced boundary config-space grid with 1 cell in dimension of self._dir.
   if self.cdim == self.grid:ndim() then
      self.confBoundaryGrid = self.boundaryGrid
   else
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      for d = 1, self.cdim do
         if d==self.bcDir then
            table.insert(reducedLower, -self.grid:dx(d)/2)
            table.insert(reducedUpper,  self.grid:dx(d)/2)
            table.insert(reducedNumCells, 1)
            table.insert(reducedCuts, 1)
         else
            table.insert(reducedLower,    self.grid:lower(d))
            table.insert(reducedUpper,    self.grid:upper(d))
            table.insert(reducedNumCells, self.grid:numCells(d))
            table.insert(reducedCuts,     self.grid:cuts(d))
         end
      end
      local reducedDecomp = CartDecomp.CartProd {
         comm      = self._splitComm,  cuts      = reducedCuts,
         writeRank = self.writeRank,   useShared = self.grid:isShared(),
      }
      self.confBoundaryGrid = Grid.RectCart {
         lower = reducedLower,  cells         = reducedNumCells,
         upper = reducedUpper,  decomposition = reducedDecomp,
      }
   end
end
function BCsBase:createDiagnostics(mySpecies) end
function BCsBase:storeBoundaryFlux(tCurr, rkIdx, qOut) end
function BCsBase:copyBoundaryFluxField(inIdx, outIdx) end
function BCsBase:combineBoundaryFluxField(outIdx, a, aIdx, ...) end
function BCsBase:computeBoundaryFluxRate(dtIn) end
function BCsBase:advance(tCurr, fIn, species, fRhsOut) end
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

return BCsBase
