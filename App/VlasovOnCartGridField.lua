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

local EmField = {}
function EmField:new(tbl)
   local self = setmetatable({}, EmField)
   self.type = "emField" -- to identify objects of this type

   self.epsilon0 = tbl.epsilon0
   self.mu0 = tbl.mu0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   -- store initial condition function
   self.initFunc = tbl.init

   self.lightSpeed = 1/math.sqrt(self.epsilon0*self.mu0)

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(EmField, { __call = function (self, o) return self.new(self, o) end })

-- methods for EM field object
EmField.__index = {
   hasEB = function (self)
      return true, false
   end,
   setCfl = function(self, cfl)
      self.cfl = cfl
   end,
   setIoMethod = function(self, ioMethod)
      self.ioMethod = ioMethod
   end,
   setBasis = function (self, basis)
      self.basis = basis
   end,
   setGrid = function(self, grid)
      self.grid = grid
   end,
   alloc = function(self)
      -- allocate fields needed in RK update
      self.em = {}
      for i = 1, 3 do
	 self.em[i] = DataStruct.Field {
	    onGrid = self.grid,
	    numComponents = self.basis:numBasis(),
	    ghost = {1, 1}
	 }
      end
      
      -- create Adios object for field I/O
      self.fieldIo = AdiosCartFieldIo {
	 elemType = self.em[1]:elemType(),
	 method = self.ioMethod,
      }
   end,
   initField = function(self)
      local project = Updater.ProjectOnBasis {
	 onGrid = self.grid,
	 basis = self.basis,
	 evaluate = self.initFunc
      }
      project:advance(0.0, 0.0, {}, {self.em[1]})
   end,
   write = function(self, frame, tm)
      if self.evolve then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
      else
	 -- if not evolving species, don't write anything except initial conditions
	 if frame == 0 then
	    self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
	 end
      end
   end,
   rkStepperFields = function(self)
      return self.em[1], self.em[2], self.em[3]
   end,
   forwardEuler = function(self, tCurr, dt, emIn, emOut)
      if self.evolve then
	 assert(false, "EmField:forwardEuler. NYI!")
      else
	 emOut:copy(emIn) -- just copy stuff over
	 return true, GKYL_MAX_DOUBLE
      end
   end,
   applyBc = function(self, tCurr, dt, emIn)
      emIn:sync()
   end,
   totalSolverTime = function(self)
      return 0.0
   end,
   volTime = function(self)
      return 0.0
   end,
   surfTime = function(self)
      return 0.0
   end,   
}

-- NoField ---------------------------------------------------------------------
--
-- Represents no field (nothing is evolved or stored)
--------------------------------------------------------------------------------

local NoField = {}
function NoField:new(tbl)
   local self = setmetatable({}, NoField)
   self.type = "noField" -- to identify objects of this type
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(NoField, { __call = function (self, o) return self.new(self, o) end })

-- methods for no field object
NoField.__index = {
   hasEB = function (self)
      return false, false
   end,
   setCfl = function(self, cfl)
   end,
   setIoMethod = function(self, ioMethod)
   end,
   setBasis = function (self, basis)
   end,
   setGrid = function(self, grid)
   end,
   alloc = function(self)
   end,
   initField = function(self)
   end,
   write = function(self, frame, tm)
   end,
   rkStepperFields = function(self)
      return nil, nil, nil
   end,
   forwardEuler = function(self, tCurr, dt, emIn, emOut)
   end,
   applyBc = function(self, tCurr, dt, emIn)
   end,
   totalSolverTime = function(self)
      return 0.0
   end,
   volTime = function(self)
      return 0.0
   end,
   surfTime = function(self)
      return 0.0
   end,   
}

return {
   EmField = EmField,
   NoField = NoField,
}
