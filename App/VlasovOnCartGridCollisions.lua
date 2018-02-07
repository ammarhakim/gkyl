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

function BgkCollisions:forwardEuler(tCurr, dt, emIn, momIn, emOut)
   local mys, mydt = self.fieldSlvr:advance(tCurr, dt, {emIn}, {emOut})
   return mys, mydt
end

   
function BgkCollisions:totalSolverTime()
   return self.fieldSlvr.totalTime + self.tmCurrentAccum
end

return {
   BgkCollisions = BgkCollisions,
}

function forwardEuler(tCurr, dt, inIdx, outIdx, speciesList)
   local spRkFields, spMomFields = {},  {}

   for i, nm in ipairs(species) do
      spOutFields[i] = speciesList[nm]:rkStepperFields()[outIdx]
      spMomFields[i] = speciesList[nm]:fluidMoments()
   end
   return updater:advance(tCurr, dt, spMomFields, spOutFields)
end
