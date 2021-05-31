-- Gkyl ------------------------------------------------------------------------
--
-- For accessing BC objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiabaticBasic          = require "App.BCs.AdiabaticBasic"
local BCsBase                 = require "App.BCs.BCsBase"
local BronoldFehskeReflection = require "App.BCs.BronoldFehskeReflection"
local GyrofluidBasic          = require "App.BCs.GyrofluidBasic"
local IncompEulerBasic        = require "App.BCs.IncompEulerBasic"
local VlasovBasic             = require "App.BCs.VlasovBasic"

return {
   AdiabaticBasic          = AdiabaticBasic,
   BCsBase                 = BCsBase,
   BronoldFehskeReflection = BronoldFehskeReflection,
   GyrofluidBasic          = GyrofluidBasic,
   IncompEulerBasic        = IncompEulerBasic,
   VlasovBasic             = VlasovBasic,
}
