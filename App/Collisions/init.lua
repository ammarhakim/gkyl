-- Gkyl ------------------------------------------------------------------------
--
-- For accessing collisions objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase    = require "App.Collisions.CollisionsBase"
local BgkCollisions     = require "App.Collisions.BgkCollisions"
local FluidDiffusion    = require "App.Collisions.FluidDiffusion"
local GkLBOCollisions   = require "App.Collisions.GkLBOCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VoronovIonization = require "App.Collisions.VoronovIonization"

return {
  CollisionsBase    = CollisionsBase,
  BgkCollisions     = BgkCollisions,
  FluidDiffusion    = FluidDiffusion,
  GkLBOCollisions   = GkLBOCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VoronovIonization = VoronovIonization,
}
