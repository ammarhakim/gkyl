-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Fluid diffusion term.
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
local Lin              = require "Lib.Linalg"
local Mpi              = require "Comm.Mpi"

-- FluidDiffusion ---------------------------------------------------------------
--
-- Add a diffusion term to the right side of a fluid equation.
--------------------------------------------------------------------------------

local FluidDiffusion = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidDiffusion:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function FluidDiffusion:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.cfl = 0.0    -- Will be replaced.

   self.diffCoeff = assert(tbl.coefficient, 
      "App.FluidDiffusion: Must specify the diffusion coefficient (vector) in 'coefficient'.")

   self.usePositivity = speciesTbl.applyPositivity    -- Use positivity preserving algorithms.

   self.tmEvalMom = 0.0
end

function FluidDiffusion:setName(nm)
   self.name = nm
end
function FluidDiffusion:setSpeciesName(nm)
   self.speciesName = nm
end

function FluidDiffusion:setCfl(cfl)
   self.cfl = cfl
end
function FluidDiffusion:setConfBasis(basis)
   self.confBasis = basis
end
function FluidDiffusion:setConfGrid(grid)
   self.confGrid = grid
end

function FluidDiffusion:createSolver()
   self.cDim      = self.confBasis:ndim()
   self.cNumBasis = self.confBasis:numBasis()

   -- Intemediate storage for output of collisions.
   self.collOut = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Zero-flux BCs.
   local zfd = { }
   for d = 1, self.cDim do
      zfd[d] = d
   end

   -- Diffusion equation.
   local constDiffusionCalc = ConstDiffusionEq {
      Dcoeff = self.diffCoeff,
      basis  = self.confBasis,
      positivity = self.usePositivity,
   }
   self.diffusionSlvr = Updater.HyperDisCont {
      onGrid             = self.confGrid,
      basis              = self.confBasis,
      cfl                = self.cfl,
      equation           = constDiffusionCalc,
--      updateDirections   = zfd, -- only update velocity directions
--      zeroFluxDirections = zfd,
   }
end

function FluidDiffusion:advance(tCurr, fIn, species, fRhsOut)

   -- Compute increment from collisions and accumulate it into output.
   self.diffusionSlvr:advance(tCurr, {fIn}, {self.collOut})

   fRhsOut:accumulate(1.0, self.collOut)

end

function FluidDiffusion:write(tm, frame)
-- Since this doesn't seem to be as big a problem in Vm as in Gk, we comment this out for now.
--   self.primMomLimitCrossings:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame)
end

function FluidDiffusion:totalTime()
   return self.diffusionSlvr.totalTime + self.tmEvalMom
end

function FluidDiffusion:slvrTime()
   return self.diffusionSlvr.totalTime
end

function FluidDiffusion:momTime()
   return self.tmEvalMom
end

return FluidDiffusion
