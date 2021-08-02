local Proto = require "Lib.Proto"

-- Empty shell collisions base class.
local CollisionsBase = Proto()

function CollisionsBase:init(tbl) self.tbl = tbl end
function CollisionsBase:fullInit(speciesTbl) end
function CollisionsBase:createDiagnostics(self, field) return nil end
function CollisionsBase:setName(nm) self.name = self.speciesName.."_"..nm end

return CollisionsBase
