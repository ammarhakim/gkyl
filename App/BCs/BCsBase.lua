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
function BCsBase:createSolver(thisSpecies, extField) end
function BCsBase:advance(tCurr, fIn, species, fRhsOut) end
function BCsBase:createDiagnostics(thisSpecies, momTable) end
function BCsBase:writeDiagnosticIntegratedMoments(tm, frame) end
function BCsBase:write(tm, frame) end

return BCsBase
