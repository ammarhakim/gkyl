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
local GkBasic                 = require "App.BCs.GkBasic"
local GyrofluidBasic          = require "App.BCs.GyrofluidBasic"
local IncompEulerBasic        = require "App.BCs.IncompEulerBasic"
local NeutralRecycling        = require "App.BCs.NeutralRecycling"
local TwistShift              = require "App.BCs.TwistShift"
local VlasovBasic             = require "App.BCs.VlasovBasic"

return {
   AdiabaticBasic          = AdiabaticBasic,
   BCsBase                 = BCsBase,
   BronoldFehskeReflection = BronoldFehskeReflection,
   GkBasic                 = GkBasic,
   GyrofluidBasic          = GyrofluidBasic,
   IncompEulerBasic        = IncompEulerBasic,
   NeutralRecyclingBasic   = NeutralRecyclingBasic,
   TwistShift              = TwistShift,
   VlasovBasic             = VlasovBasic,
}
