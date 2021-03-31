local Proto = require "Lib.Proto"

-- empty shell source base class
local SourceBase = Proto()

-- optional output for source
function SourceBase:write(tm) end

return SourceBase
