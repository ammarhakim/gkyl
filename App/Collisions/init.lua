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
local GkChargeExchange  = require "App.Collisions.GkChargeExchange"
local GkIonization      = require "App.Collisions.GkIonization"
local VmBGKCollisions   = require "App.Collisions.VmBGKCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VmChargeExchange  = require "App.Collisions.VmChargeExchange"
local VmIonization      = require "App.Collisions.VmIonization"

return {
  CollisionsBase    = CollisionsBase,
  Diffusion         = Diffusion,
  GkBGKCollisions   = GkBGKCollisions,
  GkLBOCollisions   = GkLBOCollisions,
  GkChargeExchange  = GkChargeExchange,
  GkIonization      = GkIonization,
  VmBGKCollisions   = VmBGKCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VmChargeExchange  = VmChargeExchange,
  VmIonization      = VmIonization,
}
