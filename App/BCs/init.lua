-- Gkyl ------------------------------------------------------------------------
--
-- For accessing BC objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase        = require "App.BCs.BCsBase"
local GyrofluidBasic = require "App.BCs.GyrofluidBasic"

return {
   BCsBase        = BCsBase,
   GyrofluidBasic = GyrofluidBasic,
}
