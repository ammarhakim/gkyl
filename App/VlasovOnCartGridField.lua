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
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local PerfMaxwell = require "Eq.PerfMaxwell"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "Lib.date"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end   
end

-- EmField ---------------------------------------------------------------------
--
-- Electromagnetic field (Maxwell equation are evolved)
--------------------------------------------------------------------------------

local EmField = Proto()

function EmField:init(tbl)
   self.epsilon0 = tbl.epsilon0
   self.mu0 = tbl.mu0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   self.ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self.cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0   

   -- store initial condition function (this is a wrapper around user
   -- supplied function as we need to add correction potential ICs
   -- here)
   self.initFunc = function (t, xn)
      local ex, ey, ez, bx, by, bz = tbl.init(t, xn)
      return ex, ey, ez, bx, by, bz, 0.0, 0.0
   end

   self.lightSpeed = 1/math.sqrt(self.epsilon0*self.mu0)
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

function EmField:initField()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc
   }
   project:advance(0.0, 0.0, {}, {self.em[1]})
   self:applyBc(0.0, 0.0, self.em[1])
end

function EmField:write(frame, tm)
   if self.evolve then
      self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
   else
      -- if not evolving species, don't write anything except initial conditions
      if frame == 0 then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
      end
   end
end

function EmField:rkStepperFields()
   return self.em
end

function EmField:forwardEuler(tCurr, dt, emIn, emOut)
   if self.evolve then
      return self.fieldSlvr:advance(tCurr, dt, {emIn}, {emOut})
   else
      emOut:copy(emIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function EmField:applyBc(tCurr, dt, emIn)
   emIn:sync()
end
   
function EmField:totalSolverTime()
   return self.fieldSlvr.totalTime
end

-- NoField ---------------------------------------------------------------------
--
-- Represents no field (nothing is evolved or stored)
--------------------------------------------------------------------------------

local NoField = Proto()

-- methods for no field object
function NoField:init(tbl) end
function NoField:hasEB() return false, false end
function NoField:setCfl() end
function NoField:setIoMethod(ioMethod) end
function NoField:setBasis(basis) end
function NoField:setGrid(grid) end
function NoField:alloc(nField) end
function NoField:createSolver() end
function NoField:initField() end
function NoField:write(frame, tm) end
function NoField:rkStepperFields() return {} end
function NoField:forwardEuler(tCurr, dt, emIn, emOut) return true, GKYL_MAX_DOUBLE end
function NoField:applyBc(tCurr, dt, emIn) end
function NoField:totalSolverTime() return 0.0 end

return {
   EmField = EmField,
   NoField = NoField,
}
