-- for accessing any species object
local SpeciesBase = require "PlasmaApp.Species.SpeciesBase"
local KineticSpecies = require "PlasmaApp.Species.KineticSpecies"
local VlasovSpecies = require "PlasmaApp.Species.VlasovSpecies"
--local GkSpecies = require "PlasmaApp.Species.GkSpecies"

return {
  SpeciesBase = SpeciesBase,
  KineticSpecies = KineticSpecies,
  VlasovSpecies = VlasovSpecies,
  --GkSpecies = GkSpecies,
}
