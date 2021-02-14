-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: a (hyper)diffusion term.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase   = require "App.Collisions.CollisionsBase"
local DataStruct       = require "DataStruct"
local Proto            = require "Lib.Proto"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local ConstDiffusionEq = require "Eq.ConstDiffusion"
local xsys             = require "xsys"
local Mpi              = require "Comm.Mpi"

-- Diffusion ---------------------------------------------------------------
--
-- Add a diffusion term to the right side of an equation.
--------------------------------------------------------------------------------

local Diffusion = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function Diffusion:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function Diffusion:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.diffCoeff = assert(tbl.coefficient,   -- Diffusion coefficient (or vector).
      "App.Diffusion: Must specify the diffusion coefficient (vector) in 'coefficient'.")
   self.diffDirs  = tbl.diffusiveDirs         -- Directions in which to apply diffusion.
   self.diffOrder = tbl.order                 -- Read the diffusion operator order.

   self.usePositivity = speciesTbl.applyPositivity   -- Use positivity preserving algorithms.

   self.cfl = 0.0   -- Will be replaced.

   -- Set these values to be consistent with other collision apps.
   self.collidingSpecies = self.name 
   self.selfCollisions   = false
   self.crossCollisions  = false
   self.varNu            = false
   self.timeDepNu        = false
   self.collFreqs        = {1}

   self.timers = {nonSlvr = 0.}
end

function Diffusion:setName(nm)
   self.name = nm
end
function Diffusion:setSpeciesName(nm)
   self.speciesName = nm
end
function Diffusion:setCfl(cfl)
   self.cfl = cfl
end
function Diffusion:setConfBasis(basis)
   self.confBasis = basis
end
function Diffusion:setConfGrid(grid)
   self.confGrid = grid
end
function Diffusion:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function Diffusion:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function Diffusion:createSolver()
   local grid, basis
   local zfd = {}
   if self.phaseGrid then
      -- Running a phase-space simulation.
      grid  = self.phaseGrid
      basis = self.phaseBasis
      local vdim = self.phaseGrid:ndim()-self.confGrid:ndim()
      -- Zero-flux BCs.
      for d = 1, vdim do
         zfd[d] = self.confGrid:ndim() + d
      end
   else
      -- Running a conf-space simulation.
      grid  = self.confGrid
      basis = self.confBasis
   end

   local dim = basis:ndim()
   if self.diffDirs then
      assert(#self.diffDirs<=dim, "App.Diffusion: 'diffusiveDirs' cannot have more entries than the simulation's dimensions.")
   else
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, dim do self.diffDirs[d] = d end
   end

   -- Intemediate storage for output of collisions.
   self.collOut = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   -- Diffusion equation.
   self.equation = ConstDiffusionEq {
      basis         = basis,
      coefficient   = self.diffCoeff,
      diffusiveDirs = self.diffDirs,
      positivity    = self.usePositivity,
      order         = self.diffOrder,
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = grid,
      basis              = basis,
      cfl                = self.cfl,
      equation           = self.equation,
      updateDirections   = self.diffDirs,
      zeroFluxDirections = zfd,
   }
end

function Diffusion:advance(tCurr, fIn, species, fRhsOut)

   -- Compute increment from diffusion and accumulate it into output.
   self.collisionSlvr:advance(tCurr, {fIn}, {self.collOut})

   local tmNonSlvrStart = Time.clock()
   -- Barrier over shared communicator before accumulate
   if self.phaseGrid then
      Mpi.Barrier(self.phaseGrid:commSet().sharedComm)
   end

   fRhsOut:accumulate(1.0, self.collOut)
   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end

function Diffusion:write(tm, frame)
-- Since this doesn't seem to be as big a problem in Vm as in Gk, we comment this out for now.
--   self.primMomLimitCrossings:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame)
end

function Diffusion:totalTime()
   return self.collisionSlvr.totalTime + self.timers.nonSlvr
end

function Diffusion:slvrTime()
   return self.collisionSlvr.totalTime
end

function Diffusion:nonSlvrTime()
   return self.timersNonSlvr
end

return Diffusion
