-- Gkyl ------------------------------------------------------------------------
--
-- Compute the heat flux terms in the energy and perpendicular pressure
-- equations in the gyrofluid model.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local CollisionsBase = require "App.Collisions.CollisionsBase"
local GfHeatFluxEq   = require "Eq.GyrofluidHeatFlux"
local Updater        = require "Updater"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"

local GfHeatFlux = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GfHeatFlux:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GfHeatFlux:fullInit(inputSpeciesTbl)
   local tbl = self.tbl

   self.kappaPar  = assert(tbl.kappaPar, "App.GfHeatFlux: Must specify the parallel heat diffusivity in 'kappaPar'.")
   self.kappaPerp = assert(tbl.kappaPerp, "App.GfHeatFlux: Must specify the perpendicular heat diffusivity in 'kappaPerp'.")

   self.mass, self.charge = inputSpeciesTbl.mass, inputSpeciesTbl.charge

   self.timers = {totalTime = 0.}
end

function GfHeatFlux:setName(nm) self.name = nm end
function GfHeatFlux:setSpeciesName(nm) self.speciesName = nm end
function GfHeatFlux:setConfBasis(basis) self.basis = basis end
function GfHeatFlux:setConfGrid(grid) self.grid = grid end
function GfHeatFlux:setCfl(cfl) self.cfl = cfl end

function GfHeatFlux:createSolver(mySpecies, externalField)
   local grid, basis = self.grid, self.basis
   local hfDirs = {basis:ndim()}   -- Update and boundary-surf directions (will apply).

   -- Equation object with heat flux terms.
   self.equation = GfHeatFluxEq {
      onGrid   = grid,           kappaPerp = self.kappaPerp,
      basis    = basis,          charge    = self.charge,
      kappaPar = self.kappaPar,  mass      = self.mass,      
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid = grid,      equation         = self.equation,
      basis  = basis,     updateDirections = hfDirs,
      cfl    = self.cfl,  globalUpwind     = false,
--      zeroFluxDirections = hfDirs,
   }

   -- Intemediate storage for output of this app.
   self.momDotOut = mySpecies:allocVectorMoment(mySpecies.nMoments)
end

function GfHeatFlux:advance(tCurr, momIn, species, emIn, momRhsOut)
   local tm = Time.clock()

   local em     = emIn[1]   -- Dynamic fields (e.g. phi, Apar)
   local extGeo = emIn[2]   -- Geometry/external field.

   local selfSpecies = species[self.speciesName]

   self.collisionSlvr:advance(tCurr, {momIn, em, extGeo, selfSpecies.primMomSelf}, {self.momDotOut})

   momRhsOut:accumulate(1.0, self.momDotOut)

   self.timers.totalTime = self.timers.totalTime + Time.clock() - tm
end

function GfHeatFlux:write(tm, frame, species) end

function GfHeatFlux:totalTime() return self.timers.totalTime end

function GfHeatFlux:slvrTime() return self.timers.totalTime end

function GfHeatFlux:nonSlvrTime() return 0. end

return GfHeatFlux
