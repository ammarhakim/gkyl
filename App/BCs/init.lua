-- Gkyl ------------------------------------------------------------------------
--
-- For accessing BC objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiabaticBasic          = require "App.BCs.AdiabaticBasic"
local BCsBase                 = require "App.BCs.BCsBase"
local BCsBaseGK               = require "App.BCs.BCsBaseGK"
local BCsBaseVlasov           = require "App.BCs.BCsBaseVlasov"
local BronoldFehskeReflection = require "App.BCs.BronoldFehskeReflection"
local GkBasic                 = require "App.BCs.GkBasic"
local GyrofluidBasic          = require "App.BCs.GyrofluidBasic"
local FluidBasic              = require "App.BCs.FluidBasic"
local NeutralRecycling        = require "App.BCs.NeutralRecycling"
local TwistShift              = require "App.BCs.TwistShift"
local VlasovBasic             = require "App.BCs.VlasovBasic"
local MaxwellianBc            = require "App.BCs.MaxwellianBC"

return {
   AdiabaticBasic          = AdiabaticBasic,
   BCsBase                 = BCsBase,
   BCsBaseGK               = BCsBaseGK,
   BCsBaseVlasov           = BCsBaseVlasov,
   BronoldFehskeReflection = BronoldFehskeReflection,
   GkBasic                 = GkBasic,
   GyrofluidBasic          = GyrofluidBasic,
   FluidBasic              = FluidBasic,
   NeutralRecyclingBasic   = NeutralRecyclingBasic,
   TwistShift              = TwistShift,
   VlasovBasic             = VlasovBasic,
   MaxwellianBC            = MaxwellianBC,
}
