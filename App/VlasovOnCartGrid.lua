-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Maxwell, Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local date = require "Lib.date"


local function showTbl(nm, tbl)
   io.write(string.format("--- %s ---\n", nm))
   for _, v in ipairs(tbl) do
      io.write(v); io.write(" ")
   end
   io.write("\n---\n")
end

-- For returning module table
local M = {}

-- Class to store species-specific info
local Species = {}

function Species:new(tbl)
   local self = setmetatable({}, Species)
   self.type = "species"
   self.name = "name"
   
   self.charge, self.mass = tbl.charge, tbl.mass
   self.qbym = self.charge/self.mass
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.decompCuts = tbl.decompCuts

   self.vdim = #self.cells -- velocity dimensions

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")
   assert(#self.decompCuts == self.vdim, "'decompCuts' must have " .. self.vdim .. " entries")
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Species, { __call = function (self, o) return self.new(self, o) end })

-- methods for species object
Species.__index = {
   createPhaseGrid = function(self, cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
      self.cdim = #cCells
      self.ndim = self.cdim+self.vdim

      -- create decomposition      
      local decompCuts = {}
      for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
      for d = 1, self.vdim do table.insert(decompCuts, self.decompCuts[d]) end
      self.decomp = DecompRegionCalc.CartProd {
      	 cuts = decompCuts,
      	 shared = false,
      }

      -- create computational domain
      local lower, upper, cells = {}, {}, {}
      for d = 1, self.cdim do
      	 table.insert(lower, cLo[d])
      	 table.insert(upper, cUp[d])
      	 table.insert(cells, cCells[d])
      end
      for d = 1, self.vdim do
      	 table.insert(lower, self.lower[d])
      	 table.insert(upper, self.upper[d])
      	 table.insert(cells, self.cells[d])
      end
      self.grid = Grid.RectCart {
	 lower = lower,
	 upper = upper,
	 cells = cells,
	 periodicDirs = cPeriodicDirs,
	 decomposition = self.decomp,
      }
   end,
}

-- add to module table
M.Species = Species

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = tbl.logToFile and tbl.logToFile or false
   }

   log(date(false):fmt()) -- time-stamp for sim start

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s", varNm, tostring(default)))
      return default
   end

   log("Initializing VlasovOnCartGrid simulation ...")
   local tmStart = Time.clock()

   local cdim = #tbl.lower -- configuration space dimensions
   assert(cdim == #tbl.upper, "upper should have exactly " .. cdim .. " entries")
   assert(cdim == #tbl.cells, "cells should have exactly " .. cdim .. " entries")

   -- basis function name
   local basisNm = warnDefault(tbl.basis, "basis", "serendipity")
   if basisNm ~= "serendipity" and basisNm ~= "maximal-order" then
      assert(false, "Incorrect basis type " .. basisNm .. " specified")
   end

   -- polynomial order
   local polyOrder = tbl.polyOrder
   -- CFL number
   local cfl = warnDefault(tbl.cfl, "cfl", 0.1)

   -- parallel decomposition stuff
   local decompCuts = tbl.decompCuts
   if tbl.decompCuts then
      assert(cdim == #tbl.decompCuts, "decompCuts should have exactly " .. cdim .. " entries")
   else
      -- if not specified, use 1 processor
      decompCuts = { }
      for d = 1, cdim do decompCuts[d] = 1 end
   end
   local useShared = tbl.useShared and tbl.useShared or false

   -- extract periodic directions
   local periodicDirs = {}
   if tbl.periodicDirs then
      for i, d in ipairs(tbl.periodicDirs) do
	 if d<1 and d>3 then
	    assert(false, "Directions in periodicDirs table should be 1 (for X), 2 (for Y) or 3 (for Z)")
	 end
	 periodicDirs[i] = d
      end
   end

   -- read in information about each species
   local species = { }
   for nm, val in pairs(tbl) do
      if type(val) == "table" then
	 if val.type == "species" then
	    val.name = nm
	    table.insert(species, val)
	 end
      end
   end

   -- construct phase-space grid for each species
   for _, s in ipairs(species) do
      s:createPhaseGrid(tbl.lower, tbl.upper, tbl.cells, decompCuts, periodicDirs)
   end

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n", tmEnd-tmStart))

   -- return function that runs main simulation loop   
   return function(self)
      log("Starting main loop of HyperEqnOnCartGrid simulation ...")
      local tStart, tEnd, nFrame = 0, tbl.tEnd, tbl.nFrame
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or tEnd-tStart -- initial time-step
      local frame = 1
      local tFrame = (tEnd-tStart)/nFrame
      local nextIOt = tFrame
      local step = 1
      local tCurr = tStart
      local myDt = initDt

      local tmSimStart = Time.clock()
      -- main simulation loop

      local tmSimEnd = Time.clock()
      log(string.format("Main loop completed in %g sec\n", tmSimEnd-tmSimStart))
      log(date(false):fmt()) -- time-stamp for sim end
   end
end

-- VlasovOnCartGrid application object
local App = {}
-- constructor
function App:new(tbl)
   local self = setmetatable({}, App)
   self._runApplication = buildApplication(self, tbl)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(App, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
App.__index = {
   run = function (self)
      return self:_runApplication()
   end
}

-- add to table
M.App = App

return M
