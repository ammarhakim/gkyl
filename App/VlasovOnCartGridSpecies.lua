-- Gkyl ------------------------------------------------------------------------
--
-- VlasovOnCartGrid support code: Species object
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

-- Class to store species-specific info
local Species = {}
function Species:new(tbl)
   local self = setmetatable({}, Species)
   self.type = "species" -- to identify objects of this (Species) type

   self.name = "name"
   self.cfl =  0.1
   self.charge, self.mass = tbl.charge, tbl.mass
   self.qbym = self.charge/self.mass
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim = #self.cells -- velocity dimensions
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.confBasis = nil -- Will be set later

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")

   self.decompCuts = {}
   -- parallel decomposition stuff
   if tbl.decompCuts then
      assert(self.vdim == #tbl.decompCuts, "decompCuts should have exactly " .. self.vdim .. " entries")
      self.decompCuts = tbl.decompCuts
   else
      -- if not specified, use 1 processor
      for d = 1, self.vdim do self.decompCuts[d] = 1 end
   end

   -- store initial condition function
   self.initFunc = tbl.init
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Species, { __call = function (self, o) return self.new(self, o) end })

-- methods for species object
Species.__index = {
   ndim = function(self)
      return self.vdim
   end,
   setName = function(self, nm)
      self.name = nm
   end,
   setCfl = function(self, cfl)
      self.cfl = cfl
   end,
   setIoMethod = function(self, ioMethod)
      self.ioMethod = ioMethod
   end,
   setConfBasis = function (self, basis)
      self.confBasis = basis
   end,   
   createGrid = function(self, cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
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
   createBasis = function(self, nm, polyOrder)
      self.basis = createBasis(nm, self.ndim, polyOrder)
   end,
   alloc = function(self)
      -- allocate fields needed in RK update
      self.distf = {}
      for i = 1, 3 do
	 self.distf[i] = DataStruct.Field {
	    onGrid = self.grid,
	    numComponents = self.basis:numBasis(),
	    ghost = {1, 1}
	 }
      end
      -- create Adios object for field I/O
      self.distIo = AdiosCartFieldIo {
	 elemType = self.distf[1]:elemType(),
	 method = self.ioMethod,
      }
   end,
   createSolver = function (self, hasE, hasB)
      -- create updater to advance solution by one time-step
      self.vlasovSlvr = Updater.VlasovDisCont {
	 onGrid = self.grid,
	 phaseBasis = self.basis,
	 confBasis = self.confBasis,
	 charge = self.charge,
	 mass = self.mass,
	 cfl = self.cfl,
	 hasElectricField = hasE,
	 hasMagneticField = hasB,
      }
   end,
   initDist = function(self)
      local project = Updater.ProjectOnBasis {
	 onGrid = self.grid,
	 basis = self.basis,
	 evaluate = self.initFunc
      }
      project:advance(0.0, 0.0, {}, {self.distf[1]})
   end,
   write = function(self, frame, tm)
      if self.evolve then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, frame), tm)
      else
	 -- if not evolving species, don't write anything except initial conditions
	 if frame == 0 then
	    self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, frame), tm)
	 end
      end
   end,
   rkStepperFields = function(self)
      return self.distf[1], self.distf[2], self.distf[3]
   end,
   forwardEuler = function(self, tCurr, dt, fIn, fOut)
      if self.evolve then
	 return self.vlasovSlvr:advance(tCurr, dt, {fIn}, {fOut})
      else
	 fOut:copy(fIn) -- just copy stuff over
	 return true, GKYL_MAX_DOUBLE
      end
   end,
   applyBc = function(self, tCurr, dt, fIn)
      fIn:sync()
   end,
   totalSolverTime = function(self)
      return self.vlasovSlvr.totalTime
   end,
   streamTime = function(self)
      return self.vlasovSlvr:streamTime()
   end,
   forceTime = function(self)
      return self.vlasovSlvr:forceTime()
   end,
   incrementTime = function(self)
      return self.vlasovSlvr:incrementTime()
   end
}

return Species
