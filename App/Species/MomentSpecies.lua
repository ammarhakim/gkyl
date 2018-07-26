-- Gkyl ------------------------------------------------------------------------
--
-- Species object constructed from moment equations
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local SpeciesBase = require "App.Species.SpeciesBase"
local Time = require "Lib.Time"
local Updater = require "Updater"
local xsys = require "xsys"

-- Species object treated as moment equations
local MomentSpecies = Proto(SpeciesBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function MomentSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function MomentSpecies:fullInit(appTbl)
   local tbl = self.tbl -- previously stored table

   self.name = "name"
   self.cfl =  0.1
   self.grid = nil -- will be set later
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0

   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species

   self.equation = tbl.equation -- equation system to evolve

   -- create triggers to write solutions
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- frame number for solutions

   assert(tbl.init, "Must specify initial conditions with 'init' function")
   self.initFunc = tbl.init -- pointer to initial condition functiona

   self.boundaryConditions = { } -- list of Bcs to apply
end

function MomentSpecies:getCharge() return self.charge end
function MomentSpecies:getMass() return self.mass end

function MomentSpecies:setName(nm)
   self.name = nm
end

function MomentSpecies:setCfl(cfl)
   self.cfl = cfl
end

function MomentSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end

function MomentSpecies:setConfGrid(grid)
   self.grid = grid
end

function MomentSpecies:alloc(num)
   local meqn = self.equation:numEquations()
   self.fluids = {}

   for i = 1, num do
      self.fluids[i] = DataStruct.Field {
	 onGrid = self.grid,
	 numComponents = meqn,
	 ghost = {2, 2},
      }
   end   
end

function MomentSpecies:initDist()
   local fld = self.fluids[1]
   local grid = self.grid
   local meqn = self.equation:numEquations()
   local xn = Lin.Vec(fld:ndim())
   
   local fItr = fld:get(1)      
   local indexer = fld:genIndexer()
   for idx in fld:localExtRangeIter() do
      grid:setIndex(idx)
      grid:cellCenter(xn)
      -- evaluate supplied IC function
      local v = { self.initFunc(0.0, xn) } -- braces around function put return values in table
      -- set values in cell
      fld:fill(indexer(idx), fItr) -- pointer to data in cell
      for c = 1, meqn do fItr[c] = v[c] end
   end   
end

function MomentSpecies:write(tm)
   if self.evolve then
      if self.ioTrigger(tm) then
	 self.fluids[1]:write(string.format("%s_%d.bp", self.name, self.ioFrame), tm, self.ioFrame)
	 self.ioFrame = self.ioFrame + 1
      end
   else
      if self.ioFrame == 0 then
	 self.fluids[1]:write(string.format("%s_%d.bp", self.name, self.ioFrame), tm, self.ioFrame)
	 self.ioFrame = self.ioFrame + 1
      end
   end
end

return MomentSpecies
