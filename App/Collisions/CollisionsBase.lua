local Proto = require "Lib.Proto"

-- Empty shell collisions base class.
local CollisionsBase = Proto()

function CollisionsBase:setName(nm) self.name = self.speciesName.."_"..nm end
function CollisionsBase:setCollName(nm) self.collName = nm end

return CollisionsBase
