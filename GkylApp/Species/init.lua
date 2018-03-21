-- for accessing any species object
local SpeciesBase = require "GkylApp.Species.SpeciesBase"
local KineticSpecies = require "GkylApp.Species.KineticSpecies"
local VlasovSpecies = require "GkylApp.Species.VlasovSpecies"
local GkSpecies = require "GkylApp.Species.GkSpecies"

return {
  SpeciesBase = SpeciesBase,
  KineticSpecies = KineticSpecies,
  VlasovSpecies = VlasovSpecies,
  GkSpecies = GkSpecies,
}
