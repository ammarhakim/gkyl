-- Gkyl ------------------------------------------------------------------------
--
-- For accessing species objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiabaticSpecies   = require "App.Species.AdiabaticSpecies"
local FluidSpecies       = require "App.Species.FluidSpecies"
local FuncVlasovSpecies  = require "App.Species.FuncVlasovSpecies"
local GkSpecies          = require "App.Species.GkSpecies"
local GyrofluidSpecies   = require "App.Species.GyrofluidSpecies"
local HasegawaWakataniSpecies = require "App.Species.HasegawaWakataniSpecies"
local IncompEulerSpecies = require "App.Species.IncompEulerSpecies"
local MomentSpecies      = require "App.Species.MomentSpecies"
local PassiveAdvectionSpecies = require "App.Species.PassiveAdvectionSpecies"
local SpeciesBase        = require "App.Species.SpeciesBase"
local VlasovSpecies      = require "App.Species.VlasovSpecies"

return {
   AdiabaticSpecies   = AdiabaticSpecies,
   FluidSpecies       = FluidSpecies,
   FuncVlasovSpecies  = FuncVlasovSpecies,
   GkSpecies          = GkSpecies,
   GyrofluidSpecies   = GyrofluidSpecies,
   HasegawaWakataniSpecies = HasegawaWakataniSpecies,
   IncompEulerSpecies = IncompEulerSpecies,
   MomentSpecies      = MomentSpecies,
   PassiveAdvectionSpecies = PassiveAdvectionSpecies,
   SpeciesBase        = SpeciesBase,
   VlasovSpecies      = VlasovSpecies,
}
