-- Gkyl ------------------------------------------------------------------------
--
-- Projections base object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- empty shell species base class
local ProjectionBase = Proto()

-- functions that must be defined by subclasses
function ProjectionBase:init(tbl) end
function ProjectionBase:fullInit(appTbl) end
function ProjectionBase:set() end

return ProjectionBase

