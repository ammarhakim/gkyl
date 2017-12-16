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

-- For returning module table
local M = {}

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = tbl.logToFile and tbl.logToFile or false
   }

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s", varNm, tostring(default)))
      return default
   end

   log("Initializing VlasovOnCartGrid simulation ...")
   local tmStart = Time.clock()   

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n", tmEnd-tmStart))

   -- return function that runs main simulation loop   
   return function(self)
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
