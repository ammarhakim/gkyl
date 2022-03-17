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
local TokamakEdge             = require "App.BCs.TokamakEdge"
local GyrofluidBasic          = require "App.BCs.GyrofluidBasic"
local FluidBasic              = require "App.BCs.FluidBasic"
local NeutralRecycling        = require "App.BCs.NeutralRecycling"
local TwistShift              = require "App.BCs.TwistShift"
local VlasovBasic             = require "App.BCs.VlasovBasic"

return {
   AdiabaticBasic          = AdiabaticBasic,
   BCsBase                 = BCsBase,
   BronoldFehskeReflection = BronoldFehskeReflection,
   GkBasic                 = GkBasic,
   TokamakEdge             = TokamakEdge,
   GyrofluidBasic          = GyrofluidBasic,
   FluidBasic              = FluidBasic,
   NeutralRecyclingBasic   = NeutralRecyclingBasic,
   TwistShift              = TwistShift,
   VlasovBasic             = VlasovBasic,
}
