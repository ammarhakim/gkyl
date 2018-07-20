-- Gkyl ------------------------------------------------------------------------
--
-- Species object constructed from moment equations
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local SpeciesBase = require "App.Species.SpeciesBase"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"

local MomentSpecies = Proto(SpeciesBase)

return MomentSpecies
