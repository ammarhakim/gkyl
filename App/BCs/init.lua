-- Gkyl ------------------------------------------------------------------------
--
-- For accessing BC objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiabaticBasic   = require "App.BCs.AdiabaticBasic"
local BCsBase          = require "App.BCs.BCsBase"
local GyrofluidBasic   = require "App.BCs.GyrofluidBasic"
local IncompEulerBasic = require "App.BCs.IncompEulerBasic"

return {
   AdiabaticBasic   = AdiabaticBasic,
   BCsBase          = BCsBase,
   GyrofluidBasic   = GyrofluidBasic,
   IncompEulerBasic = IncompEulerBasic,
}
