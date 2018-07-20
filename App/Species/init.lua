-- Gkyl ------------------------------------------------------------------------
--
-- For accessing species objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiabaticSpecies = require "App.Species.AdiabaticSpecies"
local FluidSpecies = require "App.Species.FluidSpecies"
local GkSpecies = require "App.Species.GkSpecies"
local IncompEulerSpecies = require "App.Species.IncompEulerSpecies"
local KineticSpecies = require "App.Species.KineticSpecies"
local MomentSpecies = require "App.Species.MomentSpecies"
local SpeciesBase = require "App.Species.SpeciesBase"
local VlasovSpecies = require "App.Species.VlasovSpecies"

return {
   AdiabaticSpecies = AdiabaticSpecies,
   FluidSpecies = FluidSpecies,
   GkSpecies = GkSpecies,
   IncompEulerSpecies = IncompEulerSpecies,
   KineticSpecies = KineticSpecies,
   MomentSpecies = MomentSpecies,
   SpeciesBase = SpeciesBase,
   VlasovSpecies = VlasovSpecies,
}
