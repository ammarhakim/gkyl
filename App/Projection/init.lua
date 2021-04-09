-- Gkyl ------------------------------------------------------------------------
--
-- For accessing projections objects.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local FluidProjection   = require "App.Projection.FluidProjection"
local GkProjection      = require "App.Projection.GkProjection"
local KineticProjection = require "App.Projection.KineticProjection"
local MomentProjection  = require "App.Projection.MomentProjection"
local ProjectionBase    = require "App.Projection.ProjectionBase"
local VlasovProjection  = require "App.Projection.VlasovProjection"

return {
   FluidProjection   = FluidProjection,
   GkProjection      = GkProjection,
   KineticProjection = KineticProjection,
   MomentProjection  = MomentProjection,
   ProjectionBase    = ProjectionBase,
   VlasovProjection  = VlasovProjection,
}
