-- Gkyl ------------------------------------------------------------------------
--
-- For accessing source objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local CollisionlessEmSource = require "App.Sources.CollisionlessEmSource"
local TenMomentRelaxSource = require "App.Sources.TenMomentRelaxSource"
local AxisymmetricMomentSource = require "App.Sources.AxisymmetricMomentSource"
local AxisymmetricPhMaxwellSource = require "App.Sources.AxisymmetricPhMaxwellSource"

return {
   SourceBase = SourceBase,
   CollisionlessEmSource = CollisionlessEmSource,
   TenMomentRelaxSource = TenMomentRelaxSource,
   AxisymmetricMomentSource = AxisymmetricMomentSource,
   AxisymmetricPhMaxwellSource = AxisymmetricPhMaxwellSource,
}
