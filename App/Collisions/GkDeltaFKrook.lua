-- Gkyl ------------------------------------------------------------------------
--
-- Krook collision operator for delta-f GK. 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local CollisionsBase = require "App.Collisions.CollisionsBase"
local Updater        = require "Updater"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"

local GkDeltaFKrook = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkDeltaFKrook:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkDeltaFKrook:fullInit(inputSpeciesTbl)
   local tbl = self.tbl

   self.collKind = "GkDeltaFKrook"    -- Type of collisions model. Useful at the species app level.

   self.nuFunc = assert(tbl.nu, "App.GkDeltaFKrook: Must specify collisionality profile with function 'nu'.")

   self.timers = {totalTime = 0}
end

function GkDeltaFKrook:setName(nm) self.name = nm end
function GkDeltaFKrook:setSpeciesName(nm) self.speciesName = nm; self.collidingSpecies = {self.speciesName} end
function GkDeltaFKrook:setConfBasis(basis) self.confBasis = basis end
function GkDeltaFKrook:setConfGrid(grid) self.confGrid = grid end
function GkDeltaFKrook:setPhaseBasis(basis) self.phaseBasis = basis end
function GkDeltaFKrook:setPhaseGrid(grid) self.phaseGrid = grid end
function GkDeltaFKrook:setCfl(cfl) self.cfl = cfl end

function GkDeltaFKrook:createSolver(externalField)

   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }
   self.nu = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }

   self.project = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = self.nuFunc, 
      onGhosts = true
   }

   self.project:advance(0., {}, {self.nu})
end

function GkDeltaFKrook:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()

   local selfSpecies = species[self.speciesName]

   selfSpecies.confPhaseWeakMultiply:advance(tCurr, {fIn, self.nu}, {self.collOut})

   fRhsOut:accumulate(-1.0, self.collOut)

   self.timers.totalTime = self.timers.totalTime + Time.clock() - tm
end

function GkDeltaFKrook:write(tm, frame, species) end

function GkDeltaFKrook:totalTime() return self.timers.totalTime end

function GkDeltaFKrook:slvrTime() return self.timers.totalTime end

function GkDeltaFKrook:nonSlvrTime() return 0. end

return GkDeltaFKrook
