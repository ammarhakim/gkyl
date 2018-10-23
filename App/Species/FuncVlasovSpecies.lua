-- Gkyl ------------------------------------------------------------------------
--
-- FuncVlasovSpecies
-- A species object with fluid fields specified as a time-dependent function
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local LinearTrigger = require "Lib.LinearTrigger"
local Proto = require "Lib.Proto"
local SpeciesBase = require "App.Species.SpeciesBase"
local Time = require "Lib.Time"
local Updater = require "Updater"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end
end

local FuncVlasovSpecies = Proto(SpeciesBase)

-- methods for no field object
function FuncVlasovSpecies:init(tbl)
   FuncVlasovSpecies.super.init(self, tbl)
   self.tbl = tbl
end

function FuncVlasovSpecies:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   self.confBasis = nil -- Will be set later

   -- create triggers to write fields
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- frame number for IO
   
   local numComponentsTable = {tbl.currentFunc(0.0, {0.0})} --dummy call to get number of components
   self.numComponents = #numComponentsTable
   -- store function to compute current drive
   self.currentFunc = function (t, xn)
      if self.numComponents == 1 then
         local Jx = tbl.currentFunc(t, xn)
         return Jx
      elseif self.numComponents == 2 then
         local Jx, Jy = tbl.currentFunc(t, xn)
         return Jx, Jy
      elseif self.numComponents == 3 then
         local Jx, Jy, Jz = tbl.currentFunc(t, xn)
         return Jx, Jy, Jz
      else
         assert(false, string.format("Number of components %d not valid", self.numComponents))
      end
   end
end

function FuncVlasovSpecies:hasEB()
   return true, true
end

function FuncVlasovSpecies:getNdim()
   return self.ndim
end
function FuncVlasovSpecies:vdim()
   return 0
end
function FuncVlasovSpecies:setName(nm)
   self.name = nm
end
function FuncVlasovSpecies:setCfl(cfl)
   self.cfl = cfl
end
function FuncVlasovSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function FuncVlasovSpecies:setConfBasis(basis)
   self.confBasis = basis
end
function FuncVlasovSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function FuncVlasovSpecies:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
   self.cdim = #cCells
   self.ndim = self.cdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = self.useShared,
   }

   -- create computational domain
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cLo[d])
      table.insert(upper, cUp[d])
      table.insert(cells, cCells[d])
   end
   self.grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cPeriodicDirs,
      decomposition = self.decomp,
   }
end

function FuncVlasovSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
end

function FuncVlasovSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis()*dim,
      ghost = {1, 1}
   }
   return m
end

function FuncVlasovSpecies:alloc(nRkDup)
   -- allocate fields needed in driven current update
   self.currentDrive = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.numComponents*self.basis:numBasis(),
      ghost = {1, 1}
   }
      
   -- create Adios object for currentDrive I/O
   self.currentDriveIo = AdiosCartFieldIo {
      elemType = self.currentDrive:elemType(),
      method = self.ioMethod,
   }   
end

function FuncVlasovSpecies:allocMomCouplingFields()
   return {currentDensity = self:allocVectorMoment(self.numComponents)}
end

function FuncVlasovSpecies:createSolver()
   self.currentSlvr = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.currentFunc
   }
end

function FuncVlasovSpecies:createDiagnostics()
end

function FuncVlasovSpecies:initField()
   self.currentSlvr:advance(0.0, 0.0, {}, {self.currentDrive})
end

-- Since calcCouplingMoments is collects the currents of all the species
-- we just need to compute the "functional" current and add it to the total current
function FuncVlasovSpecies:calcCouplingMoments(tCurr, dt, rkIdx)
   local currentOut = self:rkStepperFields()[1]
   if self.evolve then
      self.currentSlvr:advance(tCurr, dt, {}, {currentOut})
   end
end

-- current driver already has units of current
function FuncVlasovSpecies:getCharge()
   return 1.0
end

-- Maxwell's solver calls this method to get the component of a species current
function FuncVlasovSpecies:getMomDensity()
   return self.currentDrive
end

function FuncVlasovSpecies:write(tm, force)
   if self.evolve then
      if self.ioTrigger(tm) or force then
	 self.currentDriveIo:write(self.currentDrive, string.format("func_currentDrive_%d.bp", self.ioFrame), tm, self.ioFrame)
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.ioFrame == 0 then
	 self.currentDriveIo:write(self.currentDrive, string.format("func_currentDrive_%d.bp", self.ioFrame), tm, self.ioFrame)
      end
      self.ioFrame = self.ioFrame+1
   end
end

function FuncVlasovSpecies:writeRestart(tm)
   self.currentDriveIo:write(self.currentDrive, "func_currentDrive_restart.bp", tm, self.ioFrame)
end

function FuncVlasovSpecies:readRestart()
   local tm, fr = self.currentDriveIo:read(self.currentDrive, "func_currentDrive_restart.bp")
   self.currentDrive:sync() -- must get all ghost-cell data correct
   
   self.ioFrame = fr
end

function FuncVlasovSpecies:rkStepperFields()
   return { self.currentDrive, self.currentDrive, self.currentDrive, self.currentDrive }
end

function FuncVlasovSpecies:totalSolverTime()
   return self.currentSlvr.totalTime
end

return FuncVlasovSpecies
