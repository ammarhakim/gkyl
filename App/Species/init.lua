-- for accessing any species object
local SpeciesBase = require "App.Species.SpeciesBase"
local KineticSpecies = require "App.Species.KineticSpecies"
--local FluidSpecies = require "App.Species.FluidSpecies"
local VlasovSpecies = require "App.Species.VlasovSpecies"
--local GkSpecies = require "App.Species.GkSpecies"
--local IncompEulerSpecies = require "App.Species.IncompEulerSpecies"

return {
  SpeciesBase = SpeciesBase,
  KineticSpecies = KineticSpecies,
  VlasovSpecies = VlasovSpecies,
  --GkSpecies = GkSpecies,
  --IncompEulerSpecies = IncompEulerSpecies,
}
