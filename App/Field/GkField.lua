-- GkField ---------------------------------------------------------------------
--
-- App support code: Gyrokinetic fields phi and apar, solved by (perpendicular) Poisson and Ampere equations
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local DataStruct = require "DataStruct"
local LinearTrigger = require "LinearTrigger"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"
local FieldBase = require "App.Field.FieldBase"

local GkField = Proto(FieldBase.FieldBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function GkField:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function GkField:fullInit(appTbl)
   local tbl = self.tbl -- previously store table
   
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   self.isElectromagnetic = xsys.pickBool(tbl.isElectromagnetic, false) -- electrostatic by default 
   self.nPotentials = 1
   if self.isElectromagnetic then self.nPotentials = 2 end

   -- create triggers to write fields
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- frame number for IO

   -- get boundary condition settings
   -- these will be checked for consistency when the solver is initialized
   if tbl.phiBcLeft then self.phiBcLeft = tbl.phiBcLeft end
   if tbl.phiBcRight then self.phiBcRight = tbl.phiBcRight end
   if tbl.phiBcBottom then self.phiBcBottom = tbl.phiBcBottom end
   if tbl.phiBcTop then self.phiBcTop = tbl.phiBcTop end
   if tbl.aparBcLeft then self.aparBcLeft = tbl.aparBcLeft end
   if tbl.aparBcRight then self.aparBcRight = tbl.aparBcRight end
   if tbl.aparBcBottom then self.aparBcBottom = tbl.aparBcBottom end
   if tbl.aparBcTop then self.aparBcTop = tbl.aparBcTop end
   if appTbl.periodicDirs then self.periodicDirs = appTbl.periodicDirs end

   -- for storing integrated energies
   self.elecEnergy = DataStruct.DynVector { numComponents = 1 }
   self.magEnergy = DataStruct.DynVector { numComponents = 1 }

   self.tmCurrentAccum = 0.0 -- time spent in current accumulate
end

-- methods for EM field object
function GkField:hasEB() return true, self.isElectromagnetic end
function GkField:setCfl() end
function GkField:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkField:setBasis(basis) self.basis = basis end
function GkField:setGrid(grid) self.grid = grid end

function GkField:alloc(nRkDup)
   -- allocate fields needed in RK update
   -- nField is related to number of RK stages
   self.potentials = {}
   for i = 1, nRkDup do
      self.potentials[i] = {}
      self.potentials[i].phi = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
      }
      if self.isElectromagnetic then
         self.potentials[i].apar = DataStruct.Field {
               onGrid = self.grid,
               numComponents = self.basis:numBasis(),
               ghost = {1, 1}
         }
      end
   end
      
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.potentials[1].phi:elemType(),
      method = self.ioMethod,
   }

   -- create fields for total charge and current densities
   self.chargeDens = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
   }
   self.currentDens = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
   }
end

-- solve for initial fields self-consistently 
-- from initial distribution function
function GkField:initField(species)
   self:forwardEuler(0, nil, nil, species, self.potentials)
end

function GkField:rkStepperFields()
   return self.potentials
end

function GkField:createSolver()
   self.phiSlvr = Updater.FemPerpPoisson {
     onGrid = self.grid,
     basis = self.basis,
     bcLeft = self.phiBcLeft,
     bcRight = self.phiBcRight,
     bcBottom = self.phiBcBottom,
     bcTop = self.phiBcTop,
     periodicDirs = self.periodicDirs,
   }
   if self.isElectromagnetic then
     self.aparSlvr = Updater.FemPerpPoisson {
       onGrid = self.grid,
       basis = self.basis,
       bcLeft = self.aparBcLeft,
       bcRight = self.aparBcRight,
       bcBottom = self.aparBcBottom,
       bcTop = self.aparBcTop,
       periodicDirs = self.periodicDirs,
     }
   end

   self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      numComponents = 1,
      quantity = "V2"
   }

   -- need to set this flag so that field calculated self-consistently at end of full RK timestep
   self.isElliptic = true
end

function GkField:createDiagnostics()
end

function GkField:write(tm)
   if self.evolve then
      -- compute electostatic energy integrated over domain
      self.energyCalc:advance(tm, 0.0, { self.potentials[1].phi }, { self.elecEnergy })
      if self.isElectromagnetic then 
        self.energyCalc:advance(tm, 0.0, { self.potentials[1].apar }, { self.magEnergy })
      end
      
      if self.ioTrigger(tm) then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm)
         end
	 self.elecEnergy:write(string.format("elecEnergy_%d.bp", self.ioFrame), tm)
	 self.magEnergy:write(string.format("magEnergy_%d.bp", self.ioFrame), tm)
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm)
         end
      end
      self.ioFrame = self.ioFrame+1
   end
end

-- not needed for GK
function GkField:accumulateCurrent(dt, current, em)
end

function GkField:forwardEuler(tCurr, dt, potIn, species, potOut)
   if self.evolve then
      local mys2 = true
      local mydt2 = GKYL_MAX_DOUBLE
      self.chargeDens:clear(0.0)
      for nm, s in pairs(species) do
        self.chargeDens:accumulate(s:charge(), s:getDens())
      end
      local mys, mydt = self.phiSlvr:advance(tCurr, dt, {self.chargeDens}, {potOut.phi})
      if self.isElectromagnetic then
        self.currentDens:clear(0.0)
        for nm, s in pairs(species) do
          self.currentDens:accumulate(s:charge(), s:getUpar())
        end
        mys2, mydt2 = self.aparSlvr:advance(tCurr, dt, {self.currentDens}, {potOut.apar})
      end
      return mys and mys2, math.min(mydt,mydt2)
   else
      -- just copy stuff over
      potOut.phi:copy(potIn.phi)
      if self.isElectromagnetic then potOut.apar:copy(potIn.apar) end
      return true, GKYL_MAX_DOUBLE
   end
end

-- boundary conditions handled by solver. this just updates ghosts.
function GkField:applyBc(tCurr, dt, potIn)
   potIn.phi:sync()
   if self.isElectromagnetic then potIn.apar:sync() end
end
   
function GkField:totalSolverTime()
   return self.fieldSlvr.totalTime
end

function GkField:totalBcTime()
   return 0.0
end

function GkField:energyCalcTime()
   return self.energyCalc.totalTime
end

return {GkField = GkField}
