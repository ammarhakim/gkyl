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
local MultimomentConstDiffusionEq = require "Eq.MultimomentConstDiffusion"
local xsys             = require "xsys"
local Mpi              = require "Comm.Mpi"

-- Diffusion ---------------------------------------------------------------
--
-- Add a diffusion term to the right side of an equation.
--------------------------------------------------------------------------------

local Diffusion = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function Diffusion:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function Diffusion:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   -- Diffusion coefficient. This can be a single number, or a table of numbers and/or functions:
   --   1) Single number: applies that diff coefficient along all directions.
   --   2) Table of numbers: must have the same number elements as diffusiveDirs, and applies those diffusion 
   --                        coefficients along the directions indicated in diffusiveDirs.
   --   3) Functions: must be scalar functions with (t,xn) signature that specify the diffusion
   --                 coefficients along the directions indicated in diffusiveDirs. These functions
   --                 can only depend on the variable along that direction.
   self.diffCoeff = assert(tbl.coefficient,
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

function Diffusion:setName(nm)          self.name = nm end
function Diffusion:setSpeciesName(nm)   self.speciesName = nm end
function Diffusion:setCfl(cfl)          self.cfl = cfl end
function Diffusion:setConfBasis(basis)  self.confBasis = basis end
function Diffusion:setConfGrid(grid)    self.confGrid = grid end
function Diffusion:setPhaseBasis(basis) self.phaseBasis = basis end
function Diffusion:setPhaseGrid(grid)   self.phaseGrid = grid end

function Diffusion:createSolver(mySpecies, externalField)

   local nMoments = mySpecies.nMoments

   local grid, basis
   local zfd = {}
   if self.phaseGrid then
      -- Running a phase-space simulation.
      grid, basis = self.phaseGrid, self.phaseBasis
      local vdim  = self.phaseGrid:ndim()-self.confGrid:ndim()
      -- Zero-flux BCs.
      for d = 1, vdim do zfd[d] = self.confGrid:ndim() + d end
   else
      -- Running a conf-space simulation.
      grid, basis = self.confGrid, self.confBasis
   end

   local dim = basis:ndim()
   local diffCoeffType = type(self.diffCoeff)
   if self.diffDirs then
      assert(#self.diffDirs<=dim, "App.Diffusion: 'diffusiveDirs' cannot have more entries than the simulation's dimensions.")
   else
      if diffCoeffType == "table" then
         assert(#self.diffCoeff==dim, "App.Diffusion: 'coefficient' must have same number of entries as there are dimensions if not specifying 'diffusiveDirs'.")
      end
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, dim do self.diffDirs[d] = d end
   end

   self.coefficient = {}
   if diffCoeffType == "number" then
      -- Set the diffusion coefficient to the same amplitude in all directions.
      for d = 1, #self.diffDirs do self.coefficient[d] = self.diffCoeff end
   elseif diffCoeffType == "table" then
      if self.diffDirs then
         assert(#self.diffCoeff==#self.diffDirs, "App.Diffusion: 'coefficient' table must have the same number of entries as 'diffusiveDirs'.")
      else
         assert(#self.diffCoeff==dim, "App.Diffusion: 'coefficient' table must have the same number of entries as the simulation's dimensions if 'diffusiveDirs' is not given.")
      end
      self.coefficient = self.diffCoeff
   end
   -- Check if any of the diffusion coefficients are functions. If yes, project
   -- the diffusion coefficient vector onto the basis.
   local isVarCoeff = false
   for d = 1, #self.diffDirs do 
      if type(self.coefficient[d]) == "function" then isVarCoeff = true; break end
   end
   if isVarCoeff then
      -- Turn number entries into functions.
      local coeffFuncs = {}
      for d = 1, dim do coeffFuncs[d] = function(t, xn) return 0. end end
      for d = 1, #self.diffDirs do 
         if type(self.coefficient[d]) == "number" then
            coeffFuncs[self.diffDirs[d]] = function(t, zn) return self.coefficient[d] end
         else
            coeffFuncs[self.diffDirs[d]] = self.coefficient[d]
         end
      end
      local diffCoeffFunc
      if dim==1 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn) end
      elseif dim==2 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn) end
      elseif dim==3 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn), coeffFuncs[3](t,xn) end
      elseif dim==4 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn), coeffFuncs[3](t,xn),
                                               coeffFuncs[4](t,xn) end
      elseif dim==5 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn), coeffFuncs[3](t,xn),
                                               coeffFuncs[4](t,xn), coeffFuncs[5](t,xn) end
      elseif dim==6 then
         diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn), coeffFuncs[3](t,xn),
                                               coeffFuncs[4](t,xn), coeffFuncs[5](t,xn), coeffFuncs[6](t,xn) end
      end
      self.coefficient = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis()*dim,
         ghost         = {1, 1},
      }
      local projectDiffCoeff = Updater.EvalOnNodes {
         onGrid = grid,   evaluate = diffCoeffFunc,
         basis  = basis,  onGhosts = true
      }
      projectDiffCoeff:advance(0.0, {}, {self.coefficient})
   end

   -- Intemediate storage for output of collisions.
   local compV = nMoments and nMoments or 1
   self.collOut = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*compV,
      ghost         = {1, 1},
   }
   -- Diffusion equation.
   if nMoments ~= nil and nMoments > 1 then
      self.equation = MultimomentConstDiffusionEq {
         basis         = basis,             positivity = self.usePositivity,
         coefficient   = self.coefficient,  order      = self.diffOrder,
         diffusiveDirs = self.diffDirs,     numMoments = nMoments,
      }
   else
      self.equation = ConstDiffusionEq {
         basis         = basis,             positivity = self.usePositivity,
         coefficient   = self.coefficient,  order      = self.diffOrder,
         diffusiveDirs = self.diffDirs,
      }
   end
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid = grid,      equation           = self.equation,
      basis  = basis,     updateDirections   = self.diffDirs,
      cfl    = self.cfl,  zeroFluxDirections = zfd,
   }
end

function Diffusion:advance(tCurr, fIn, species, fRhsOut)

   -- Compute increment from diffusion and accumulate it into output.
   self.collisionSlvr:advance(tCurr, {fIn}, {self.collOut})

   local tmNonSlvrStart = Time.clock()
   -- Barrier over shared communicator before accumulate
   if self.phaseGrid then Mpi.Barrier(self.phaseGrid:commSet().sharedComm) end

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
   return self.timers.nonSlvr
end

return Diffusion
