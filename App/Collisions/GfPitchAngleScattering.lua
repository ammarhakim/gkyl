-- Gkyl ------------------------------------------------------------------------
--
-- Add terms to the perpendicular pressure equation in the gyrofluid model
-- to account for pitch angle scattering.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local CollisionsBase = require "App.Collisions.CollisionsBase"
local Time           = require "Lib.Time"

local GfPitchAngleScattering = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GfPitchAngleScattering:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GfPitchAngleScattering:fullInit(inputSpeciesTbl)
   local tbl = self.tbl

   self.collidingSpecies = assert(tbl.collideWith, "App.GfPitchAngleScattering: Must specify names of species to collide with in 'collideWith'.")
   self.collFreqs        = assert(tbl.frequencies, "App.GfPitchAngleScattering: Must specify the collision frequencies in 'frequencies'.")

   self.collFreqSelf = self.collFreqs[1]

   self.timers = {totalTime = 0.}
end

function GfPitchAngleScattering:setName(nm) self.name = nm end
function GfPitchAngleScattering:setSpeciesName(nm) self.speciesName = nm end
function GfPitchAngleScattering:setConfBasis(basis) self.basis = basis end
function GfPitchAngleScattering:setConfGrid(grid) self.grid = grid end
function GfPitchAngleScattering:setCfl(cfl) self.cfl = cfl end

function GfPitchAngleScattering:createSolver(mySpecies, externalField)
   self.weakMultiply = mySpecies.weakMultiply

   self.collOut = mySpecies:allocMoment()
   self.unitFld1comp = mySpecies:allocCartField(self.grid, 1, {1,1})
   self.unitFld1comp:clear(1.)

   -- Inverse of background magnetic field.
   self.bmagInv = externalField.geo.bmagInv

   -- The species looks for a solver object with a setDtAndCflRate method. Mimic it.
   -- MF 2021/04/22: I think setDtAndCflRate should disappear.
   self.collisionSlvr = {}
   function self.collisionSlvr:setDtAndCflRate(dtGlobal, cflRateByCell)
      self.cflRateByCell = cflRateByCell
   end
end

function GfPitchAngleScattering:advance(tCurr, momIn, species, emIn, momRhsOut)
   local tm = Time.clock()

   local selfSpecies       = species[self.speciesName]
   local pJac              = selfSpecies:getPressuresJac()
   local pParJac, pPerpJac = pJac[1], pJac[2]

   self.collOut:combine(1., pParJac, -1., pPerpJac)
   self.weakMultiply:advance(tCurr, {self.bmagInv, self.collOut}, {self.collOut})

   momRhsOut:accumulateOffset(self.collFreqSelf, self.collOut, selfSpecies:getMomOff(4))

   -- Set the CFL rate from this operator.
   self.collisionSlvr.cflRateByCell:accumulate(self.collFreqSelf, self.unitFld1comp)

   self.timers.totalTime = self.timers.totalTime + Time.clock() - tm
end

function GfPitchAngleScattering:write(tm, frame, species) end

function GfPitchAngleScattering:totalTime() return self.timers.totalTime end

function GfPitchAngleScattering:slvrTime() return self.timers.totalTime end

function GfPitchAngleScattering:nonSlvrTime() return 0. end

return GfPitchAngleScattering
