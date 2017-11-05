-- Gkyl ------------------------------------------------------------------------
--
-- Euler solver on a Cartesian grid. Works in 1D, 2D and 3D.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local HyperEquation = require "Eq"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Time = require "Lib.Time"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

-- top-level method to build simulation "run" method
local function buildSimulation(self, tbl)
   local log = Logger { }  -- logger

   -- basic parameters and checks
   local ndim = #tbl.lower -- simulation dimension
   assert(ndim == #tbl.upper, "upper should have exactly " .. ndim .. " entries")
   assert(ndim == #tbl.cells, "cells should have exactly " .. ndim .. " entries")

   -- parallel decomposition stuff
   local decompCuts = tbl.decompCuts
   if tbl.decompCuts then
      assert(ndim == #tbl.decompCuts, "decompCuts should have exactly " .. ndim .. " entries")
   else
      -- if not specified, use 1 processor
      decompCuts = { }
      for d = 1, ndim do decompCuts[d] = 1 end
   end
   local useShared = tbl.useShared and tbl.useShared or false

   -- create decomposition and grid
   local decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = useShared
   }
   -- computational domain
   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      decomposition = decomp,
   }

   -- solution
   fluid = DataStruct.Field {
      onGrid = grid,
      numComponents = 5,
      ghost = {2, 2},
   }

   -- create ADIO object for field I/O
   local fieldIo = AdiosCartFieldIo { fluid:elemType() }

   -- function to apply initial conditions
   local function init(field)
      local xn = Lin.Vec(ndim)
      local indexer = field:genIndexer()
      for idx in field:localRangeIter() do
	 -- get cell-center coordinates
	 grid:setIndex(idx)
	 grid:cellCenter(xn)
	 -- evaluate supplied IC function
	 local rho, rhou, rhov, rhow, E = tbl.init(0.0, xn)

	 -- set values in cell
	 local fItr = field:get(indexer(idx)) -- pointer to data in cell
	 fItr[1], fItr[2], fItr[3], fItr[4], fItr[5] = rho, rhou, rhov, rhow, E
      end
   end
   -- initialize field
   init(fluid)
   -- write out ICs (allows user to debug sim without running it)
   fieldIo:write(fluid, "fluid_0.bp", 0.0)   

   return function (self)

   end
end

-- create default logger
log = Logger { }

-- Euler simulation object
local Sim = {}
-- constructor
function Sim:new(tbl)
   local self = setmetatable({}, Sim)
   self._runSimulation = buildSimulation(self, tbl)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Sim, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
Sim.__index = {
   run = function (self)
      return self._runSimulation(self)
   end
}

return {
   Sim = Sim,
}
