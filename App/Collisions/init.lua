-- Gkyl ------------------------------------------------------------------------
--
-- For accessing collisions objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase    = require "App.Collisions.CollisionsBase"
local FluidDiffusion    = require "App.Collisions.FluidDiffusion"
local GkBGKCollisions   = require "App.Collisions.GkBGKCollisions"
local GkLBOCollisions   = require "App.Collisions.GkLBOCollisions"
local VmBGKCollisions   = require "App.Collisions.VmBGKCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VoronovIonization = require "App.Collisions.VoronovIonization"

return {
  CollisionsBase    = CollisionsBase,
  FluidDiffusion    = FluidDiffusion,
  GkBGKCollisions   = GkBGKCollisions,
  GkLBOCollisions   = GkLBOCollisions,
  VmBGKCollisions   = VmBGKCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VoronovIonization = VoronovIonization,
}
