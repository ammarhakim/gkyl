local Proto = require "Lib.Proto"

-- empty shell source base class
local FluidSourceBase = Proto()

-- optional output for source
function FluidSourceBase:write(tm) end

return FluidSourceBase
