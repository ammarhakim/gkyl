-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the equations modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local xsys = require "xsys"
local Euler = require "Eq.Euler"
local PhMaxwell = require "Eq.PerfMaxwell"
local TenMoment = require "Eq.TenMoment"
local Advection = require "Eq.Advection"
local Burgers = require "Eq.Burgers"

return xsys.table.union(Euler, PhMaxwell, TenMoment, Advection, Burgers)
