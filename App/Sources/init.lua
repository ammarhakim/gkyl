-- Gkyl ------------------------------------------------------------------------
--
-- For accessing source objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local CollisionlessEmSource = require "App.Sources.CollisionlessEmSource"
local TenMomentRelaxSource = require "App.Sources.TenMomentRelaxSource"

return {
   SourceBase = SourceBase,
   CollisionlessEmSource = CollisionlessEmSource,
   TenMomentRelaxSource = TenMomentRelaxSource
}
