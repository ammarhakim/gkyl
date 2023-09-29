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
local GkSheath                = require "App.BCs.GkSheath"
local GyrofluidBasic          = require "App.BCs.GyrofluidBasic"
local FluidBasic              = require "App.BCs.FluidBasic"
local NeutralRecycling        = require "App.BCs.NeutralRecycling"
local TwistShift              = require "App.BCs.TwistShift"
local VlasovBasic             = require "App.BCs.VlasovBasic"
local VlasovEmission          = require "App.BCs.VlasovEmission"
local GkMaxwellianBc          = require "App.BCs.GkMaxwellianBC"

return {
   AdiabaticBasic          = AdiabaticBasic,
   BCsBase                 = BCsBase,
   BronoldFehskeReflection = BronoldFehskeReflection,
   FluidBasic              = FluidBasic,
   GkBasic                 = GkBasic,
   GkSheath                = GkSheath,
   GyrofluidBasic          = GyrofluidBasic,
   FluidBasic              = FluidBasic,
   GkMaxwellianBC          = GkMaxwellianBC,
   NeutralRecyclingBasic   = NeutralRecyclingBasic,
   TwistShift              = TwistShift,
   VlasovBasic             = VlasovBasic,
   VlasovEmission          = VlasovEmission,
}
