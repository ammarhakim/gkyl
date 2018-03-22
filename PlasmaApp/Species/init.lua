-- for accessing any species object
local SpeciesBase = require "PlasmaApp.Species.SpeciesBase"
local KineticSpecies = require "PlasmaApp.Species.KineticSpecies"
local FluidSpecies = require "PlasmaApp.Species.FluidSpecies"
local VlasovSpecies = require "PlasmaApp.Species.VlasovSpecies"
--local GkSpecies = require "PlasmaApp.Species.GkSpecies"
local IncompEulerSpecies = require "PlasmaApp.Species.IncompEulerSpecies"

return {
  SpeciesBase = SpeciesBase,
  KineticSpecies = KineticSpecies,
  VlasovSpecies = VlasovSpecies,
  --GkSpecies = GkSpecies,
  IncompEulerSpecies = IncompEulerSpecies,
}
