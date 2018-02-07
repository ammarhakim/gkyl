-- Gkyl ------------------------------------------------------------------------
--
-- VlasovOnCartGrid support code: Various field objects
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local PerfMaxwell = require "Eq.PerfMaxwell"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "Lib.date"
local xsys = require "xsys"

-- BgkCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator
--------------------------------------------------------------------------------

local BgkCollisions = Proto()

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function BgkCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function BgkCollisions:fullInit(vlasovTbl)
   local tbl = self.tbl -- previously store table
end

-- methods for Bgk collisions object
function BgkCollisions:setCfl(cfl) self.cfl = cfl end
function BgkCollisions:setBasis(basis) self.basis = basis end
function BgkCollisions:setGrid(grid) self.grid = grid end

function BgkCollisions:alloc(nField)
   -- allocate fields needed in RK update
   self.em = {}
   for i = 1, nField do
      self.em[i] = DataStruct.Field {
	 onGrid = self.grid,
	 numComponents = 8*self.basis:numBasis(),
	 ghost = {1, 1}
      }
   end
      
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.em[1]:elemType(),
      method = self.ioMethod,
   }
end

function BgkCollisions:createSolver()
   local maxwellEqn = PerfMaxwell {
      lightSpeed = self.lightSpeed,
      elcErrorSpeedFactor = self.ce,
      mgnErrorSpeedFactor = self.cb,
      basis = self.basis,
   }
   
   self.fieldSlvr = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = maxwellEqn
   }

   self.emEnergyCalc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      numComponents = 8,
      quantity = "V2"
   }
end



function BgkCollisions:rkStepperFields()
   return self.em
end

function BgkCollisions:accumulateCurrent(dt, current, em)
   if current == nil then return end

   local tmStart = Time.clock()

   -- these many current components are supplied
   local cItr, eItr = current:get(1), em:get(1)
   local cIdxr, eIdxr = current:genIndexer(), em:genIndexer()

   for idx in em:localRangeIter() do
      current:fill(cIdxr(idx), cItr)
      em:fill(eIdxr(idx), eItr)
      for i = 1, current:numComponents() do
   	 eItr[i] = eItr[i]-dt/self.epsilon0*cItr[i]
      end
   end
   self.tmCurrentAccum = self.tmCurrentAccum + Time.clock()-tmStart
end

-- momIn[1] is the current density
function BgkCollisions:forwardEuler(tCurr, dt, emIn, momIn, emOut)
   if self.evolve then
      local mys, mydt = self.fieldSlvr:advance(tCurr, dt, {emIn}, {emOut})
      self:accumulateCurrent(dt, momIn[1], emOut)
      return mys, mydt
   else
      emOut:copy(emIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function BgkCollisions:applyBc(tCurr, dt, emIn)
   emIn:sync()
end
   
function BgkCollisions:totalSolverTime()
   return self.fieldSlvr.totalTime + self.tmCurrentAccum
end

function BgkCollisions:energyCalcTime()
   return self.emEnergyCalc.totalTime
end

-- NoCollisions ----------------------------------------------------------------
--
-- Represents no collisions (nothing is evolved or stored)
--------------------------------------------------------------------------------

local NoCollisions = Proto()

-- methods for no field object
function NoCollisions:init(tbl) end
function NoCollisions:fullInit(tbl) end
function NoCollisions:setCfl() end
function NoCollisions:setBasis(basis) end
function NoCollisions:setGrid(grid) end
function NoCollisions:alloc(nField) end
function NoCollisions:createSolver() end
function NoCollisions:write(tm) end
function NoCollisions:rkStepperFields() return {} end
function NoCollisions:forwardEuler(tCurr, dt, momIn, gIn, fOut) return true, GKYL_MAX_DOUBLE end
function NoCollisions:totalSolverTime() return 0.0 end

return {
   BgkCollisions = BgkCollisions,
   NoCollisions = NoCollisions,
}
