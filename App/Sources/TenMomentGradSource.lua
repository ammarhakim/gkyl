-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Collisionless EM sources for use in fluid sims
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"

-- TenMomentGradSource ---------------------------------------------------------------
--
-- Coupled Lorentz and current sources
--------------------------------------------------------------------------------

local TenMomentGradSource = Proto(SourceBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function TenMomentGradSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function TenMomentGradSource:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   -- set later by various set() methods
   self.cfl = nil
   self.grid = nil
   self.slvr = nil

   self.speciesList = tbl.species -- list of species to update
   -- factor multiplying the thermal conductivity
   self.alpha = tbl.alpha
   -- factor to determine number of stages in super-time-stepping
   self.dtRatio = tbl.dtRatio
end

function TenMomentGradSource:setName(nm)
   self.name = nm
end

function TenMomentGradSource:setCfl(cfl)
   self.cfl = cfl
end
function TenMomentGradSource:setConfGrid(grid)
   self.grid = grid
end

function TenMomentGradSource:createSolver(species, field)
   local numSpecies = #self.speciesList
   local mass, charge = {}, {}

   local source_type
   for i, nm in ipairs(self.speciesList) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
      if not source_type then
         source_type = species[nm].nMoments
      else
         -- FIXME Currently all species must have the same moments.
         assert(source_type == species[nm].nMoments)
      end
   end

   self.slvr = Updater.TenMomentGrad {
      onGrid = self.grid,
      alpha = self.alpha,
   }
   self.ConservToPrim = Updater.TenMomentConservToPrim {
      onGrid = self.grid,
      mode = "ToPrimitive",
   }
   self.PrimToConserv = Updater.TenMomentConservToPrim {
      onGrid = self.grid,
      mode = "ToConservative",
   }
   -- Array with ten components per cell to store symmetrized heat flux.
   self.q = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   self.q:clear(0.0)
   -- Arrays for super-time-stepping
   self.qDiff0 = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   self.qDiff = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   self.qJ1 = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   self.qJ2 = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   -- six component array for advancing solution in time without affecting other components
   self.qJ = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 6,
      ghost = {2, 2},
   }
   self.qDiff0:clear(0.0)
   self.qDiff:clear(0.0)
   self.qJ:clear(0.0)
   self.qJ1:clear(0.0)
   self.qJ2:clear(0.0)
end

function TenMomentGradSource:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   return status, dtSuggested
end

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594626. doi:10.1016/j.jcp.2013.08.021)
function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
--
function a(j) return 1-b(j) end
--
function w1(s) return 4/(s^2+s-2) end
-- 
function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end
--
function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
-- 
function nu(j) return -(j-1)/j*b(j)/b(j-2) end
--
function gbar(s,j) return -a(j-1)*mubar(s,j) end
-- 
function calcNumStages(dhdp) 
   return math.ceil (math.sqrt(4*dhdp+9/4) - 1/2)
end

-- update the pressure tensor to the new time step with super-time-stepping
function TenMomentGradSource:updateSource(tCurr, dt, speciesVar, fieldVar)
   local qIn
   -- Get the name of the species being updates and store that in the outVars table
   for i, nm in ipairs(self.speciesList) do
      qIn = speciesVar[nm]
   end
   self.slvr:setDtAndCflRate(dt, nil)
   -- Convert second tensor moment from full stress tensor to just pressure tensor
   self.ConservToPrim:advance(tCurr, {}, {qIn})
   local numStages = calcNumStages(self.dtRatio)
   qIn:sync()
   -- we need this in each stage
   self.slvr:advance(tCurr, {qIn, self.q}, {self.qDiff0})

   -- stage 1
   self.qJ2:copy(qIn)
   self.qJ1:combine(1.0, qIn, mubar(numStages,1)*dt, self.qDiff0)

   -- rest of stages
   for j = 2, numStages do
      self.qJ1:sync()
      self.slvr:advance(tCurr, {self.qJ1, self.q}, {self.qDiff})
      self.qJ:combineOffset(mu(j), self.qJ1, 4, nu(j), self.qJ2, 4, 1-mu(j)-nu(j), qIn, 4,
      mubar(numStages,j)*dt, self.qDiff, 4, gbar(numStages,j)*dt, self.qDiff0, 4)

      -- reset fields for next stage
      self.qJ2:copy(self.qJ1)
      self.qJ1:copyOffset(self.qJ, 4)
   end
   qIn:copy(self.qJ1)
   -- Add the Reynolds stress (n u_i u_j) back into the result
   self.PrimToConserv:advance(tCurr, {}, {qIn})   
   qIn:sync()
   return true, GKYL_MAX_DOUBLE
end

function TenMomentGradSource:write(tm, frame)
end

function TenMomentGradSource:totalTime()
   return self.slvr.totalTime
end

return TenMomentGradSource