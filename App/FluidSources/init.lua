-- Gkyl ------------------------------------------------------------------------
--
-- For accessing fluid source objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local FluidSourceBase = require "App.FluidSources.FluidSourceBase"
local CollisionlessEmSource = require "App.Sources.CollisionlessEmSource"
local TenMomentRelaxSource = require "App.Sources.TenMomentRelaxSource"
local TenMomentGradSource = require "App.Sources.TenMomentGradSource"
local AxisymmetricMomentSource = require "App.Sources.AxisymmetricMomentSource"
local AxisymmetricPhMaxwellSource = require "App.Sources.AxisymmetricPhMaxwellSource"
local BraginskiiHeatConductionSource = require "App.Sources.BraginskiiHeatConductionSource"
local BraginskiiViscosityDiffusionSource = require "App.Sources.BraginskiiViscosityDiffusionSource"

return {
   FluidSourceBase = FluidSourceBase,
   CollisionlessEmSource = CollisionlessEmSource,
   TenMomentRelaxSource = TenMomentRelaxSource,
   TenMomentGradSource = TenMomentGradSource,
   AxisymmetricMomentSource = AxisymmetricMomentSource,
   AxisymmetricPhMaxwellSource = AxisymmetricPhMaxwellSource,
   BraginskiiHeatConductionSource = BraginskiiHeatConductionSource,
   BraginskiiViscosityDiffusionSource = BraginskiiViscosityDiffusionSource,
}
