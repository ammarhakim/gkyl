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

-- EmField ---------------------------------------------------------------------
--
-- Electromagnetic field (Maxwell equation are evolved)
--------------------------------------------------------------------------------

local EmField = Proto()

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function EmField:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function EmField:fullInit(vlasovTbl)
   local tbl = self.tbl -- previously store table
   
   self.epsilon0 = tbl.epsilon0
   self.mu0 = tbl.mu0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   self.ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self.cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0

   -- create triggers to write fields
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, vlasovTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, vlasovTbl.tEnd, vlasovTbl.nFrame)
   end

   self.ioFrame = 0 -- frame number for IO

   -- store initial condition function (this is a wrapper around user
   -- supplied function as we need to add correction potential ICs
   -- here)
   self.initFunc = function (t, xn)
      local ex, ey, ez, bx, by, bz = tbl.init(t, xn)
      return ex, ey, ez, bx, by, bz, 0.0, 0.0
   end

   self.lightSpeed = 1/math.sqrt(self.epsilon0*self.mu0)

   self.tmCurrentAccum = 0.0 -- time spent in current accumulate
end

-- methods for EM field object
function EmField:hasEB() return true, false end
function EmField:setCfl(cfl) self.cfl = cfl end
function EmField:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function EmField:setBasis(basis) self.basis = basis end
function EmField:setGrid(grid) self.grid = grid end

function EmField:alloc(nField)
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

function EmField:createSolver()
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
end

function EmField:createDiagnostics()
end

function EmField:initField()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc
   }
   project:advance(0.0, 0.0, {}, {self.em[1]})
   self:applyBc(0.0, 0.0, self.em[1])
end

function EmField:write(tm)
   if self.evolve then
      if self.ioTrigger(tm) then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", self.ioFrame), tm)
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", self.ioFrame), tm)
      end
      self.ioFrame = self.ioFrame+1
   end
end

function EmField:rkStepperFields()
   return self.em
end

function EmField:accumulateCurrent(dt, current, em)
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
function EmField:forwardEuler(tCurr, dt, emIn, momIn, emOut)
   if self.evolve then
      local mys, mydt = self.fieldSlvr:advance(tCurr, dt, {emIn}, {emOut})
      self:accumulateCurrent(dt, momIn[1], emOut)
      return mys, mydt
   else
      emOut:copy(emIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function EmField:applyBc(tCurr, dt, emIn)
   emIn:sync()
end
   
function EmField:totalSolverTime()
   return self.fieldSlvr.totalTime + self.tmCurrentAccum
end

-- NoField ---------------------------------------------------------------------
--
-- Represents no field (nothing is evolved or stored)
--------------------------------------------------------------------------------

local NoField = Proto()

-- methods for no field object
function NoField:init(tbl) end
function NoField:fullInit(tbl) end
function NoField:hasEB() return false, false end
function NoField:setCfl() end
function NoField:setIoMethod(ioMethod) end
function NoField:setBasis(basis) end
function NoField:setGrid(grid) end
function NoField:alloc(nField) end
function NoField:createSolver() end
function NoField:createDiagnostics() end
function NoField:initField() end
function NoField:write(tm) end
function NoField:rkStepperFields() return {} end
function NoField:forwardEuler(tCurr, dt, momIn, emIn, emOut) return true, GKYL_MAX_DOUBLE end
function NoField:applyBc(tCurr, dt, emIn) end
function NoField:totalSolverTime() return 0.0 end

return {
   EmField = EmField,
   NoField = NoField,
}
