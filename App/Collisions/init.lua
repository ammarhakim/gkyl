-- Gkyl ------------------------------------------------------------------------
--
-- For accessing collisions objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase    = require "App.Collisions.CollisionsBase"
local Diffusion         = require "App.Collisions.Diffusion"
local GkBGKCollisions   = require "App.Collisions.GkBGKCollisions"
local GkLBOCollisions   = require "App.Collisions.GkLBOCollisions"
local VmBGKCollisions   = require "App.Collisions.VmBGKCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VoronovIonization = require "App.Collisions.VoronovIonization"

return {
  CollisionsBase    = CollisionsBase,
  Diffusion         = Diffusion,
  GkBGKCollisions   = GkBGKCollisions,
  GkLBOCollisions   = GkLBOCollisions,
  VmBGKCollisions   = VmBGKCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VoronovIonization = VoronovIonization,
}
