-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the equations modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local xsys = require "xsys"
local Euler = require "Eq.Euler"
local PhMaxwell = require "Eq.PerfMaxwell"

return xsys.table.union(Euler, PhMaxwell)
