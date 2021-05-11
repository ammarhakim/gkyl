-- Gkyl ------------------------------------------------------------------------
--
-- Boundary condition base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell source base class.
local BCsBase = Proto()

-- Functions that must be defined by subclasses.
function BCsBase:init(tbl) self.tbl = tbl end
function BCsBase:fullInit(speciesTbl) end
function BCsBase:setSpeciesName(nm) self.speciesName = nm end
function BCsBase:setConfBasis(basis) self.basis = basis end
function BCsBase:setConfGrid(grid) self.grid = grid end
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
function BCsBase:advance(tCurr, fIn, species, fRhsOut) end
function BCsBase:storeBoundaryFlux(tCurr, rkIdx, qOut) end
function BCsBase:copyBoundaryFluxField(inIdx, outIdx) end
function BCsBase:combineBoundaryFluxField(outIdx, a, aIdx, ...) end
function BCsBase:getDir() return self.bcDir end
function BCsBase:getEdge() return self.bcEdge end
function BCsBase:getKind() return self.bcKind end

return BCsBase
