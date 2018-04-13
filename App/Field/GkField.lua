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
local AdiabaticSpecies = require "App.Species.AdiabaticSpecies"

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

   self.adiabatic = false
   if tbl.adiabatic then
      self.adiabatic = true
      self.adiabDens = tbl.adiabatic.dens
      self.adiabCharge = tbl.adiabatic.charge
      self.adiabQneutFac = tbl.adiabatic.dens*tbl.adiabatic.charge^2/tbl.adiabatic.temp
   end
   self.discontinuous = xsys.pickBool(tbl.discontinuous, false)
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
      self.potentials[i].phi:clear(0.0)
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
   self:forwardEuler(0, nil, nil, species, self.potentials[1])
end

function GkField:rkStepperFields()
   return self.potentials
end

function GkField:createSolver()
   if self.adiabatic then
      assert(self.grid:ndim()==2, "Adiabatic solve only supported in 2d for now!")
      -- if adiabatic, we don't solve the Poisson equation, but
      -- we do need to project the discontinuous charge densities
      -- onto a continuous basis to set phi. 
      -- HACK for 2d: use FemPerpPoisson with no laplacian and a weighted mass matrix.
      self.phiSlvr = Updater.FemPerpPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcLeft = self.phiBcLeft,
        bcRight = self.phiBcRight,
        bcBottom = self.phiBcBottom,
        bcTop = self.phiBcTop,
        periodicDirs = self.periodicDirs,
        laplacianWeight = 0.0, 
        modifierConstant = self.adiabQneutFac -- = n0*q^2/T
      }
   else
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
	 if self.isElectromagnetic then self.magEnergy:write(string.format("magEnergy_%d.bp", self.ioFrame), tm) end
	 
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
   if self.evolve or tCurr == 0.0 then
      local mys, mys2 = true, true
      local mydt, mydt2 = GKYL_MAX_DOUBLE, GKYL_MAX_DOUBLE
      self.chargeDens:clear(0.0)
      for nm, s in pairs(species) do
         self.chargeDens:accumulate(s:getCharge(), s:getDens())
      end
      if self.adiabatic and self.discontinuous then
         potOut.phi:combine(1.0/self.adiabQneutFac, self.chargeDens)
      else
         local mys, mydt = self.phiSlvr:advance(tCurr, dt, {self.chargeDens}, {potOut.phi})
         if self.isElectromagnetic then
           self.currentDens:clear(0.0)
           for nm, s in pairs(species) do
             self.currentDens:accumulate(s:getCharge(), s:getUpar())
           end
           mys2, mydt2 = self.aparSlvr:advance(tCurr, dt, {self.currentDens}, {potOut.apar})
         end
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
   potIn.phi:sync(true)
   if self.isElectromagnetic then potIn.apar:sync() end
end
   
function GkField:totalSolverTime()
   if self.phiSlvr then
     local time = self.phiSlvr.totalTime
     if self.isElectromagnetic and self.aparSlvr then time = time + self.aparSlvr.totalTime end
     return time
   end
   return 0.0
end

function GkField:totalBcTime()
   return 0.0
end

function GkField:energyCalcTime()
   return self.energyCalc.totalTime
end

-- GkGeometry ---------------------------------------------------------------------
--
-- A field object with fields specifying the magnetic geometry for GK
--------------------------------------------------------------------------------

local GkGeometry = Proto(FieldBase.FuncFieldBase)

-- methods for no field object
function GkGeometry:init(tbl)
   self.tbl = tbl
end

function GkGeometry:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, false) -- by default these fields are not time-dependent

   -- create triggers to write fields
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- frame number for IO
   
   -- get function to initialize background magnetic field
   self.bmagFunc = assert(tbl.bmag, "GkGeometry: must specify background magnetic field with 'bmag'")
   -- get function to initialize bcurvY = 1/B*curl(bhat).grad(y)
   self.bcurvYFunc = tbl.bcurvY

   -- wall potential for sheath BCs
   self.phiWallFunc = tbl.phiWall
end

function GkGeometry:hasEB() end
function GkGeometry:setCfl() end
function GkGeometry:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkGeometry:setBasis(basis) self.basis = basis end
function GkGeometry:setGrid(grid) self.grid = grid end

function GkGeometry:alloc()
   -- allocate fields 
   self.geo = {}

   -- background magnetic field
   self.geo.bmag = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- 1/B
   self.geo.bmagInv = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- 1/B*curl(bhat).grad(y) ... contains curvature info
   self.geo.bcurvY = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- wall potential for sheath BCs
   self.geo.phiWall = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
      
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.geo.bmag:elemType(),
      method = self.ioMethod,
   }   
end

function GkGeometry:createSolver()
   self.setBmag = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.bmagFunc,
      projectOnGhosts = true,
   }
   local bmagInvFunc = function (t, xn)
      return 1/self.bmagFunc(t,xn)
   end
   self.setBmagInv = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = bmagInvFunc
   }
   if self.bcurvYFunc then 
      self.setBcurvY = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         projectOnGhosts = true,
         evaluate = self.bcurvYFunc
      }
   end
   if self.phiWallFunc then 
      self.setPhiWall = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         projectOnGhosts = true,
         evaluate = self.phiWallFunc
      }
   end
end

function GkGeometry:createDiagnostics()
end

function GkGeometry:initField()
   self.setBmag:advance(0.0, 0.0, {}, {self.geo.bmag})
   self.setBmagInv:advance(0.0, 0.0, {}, {self.geo.bmagInv})
   if self.setBcurvY then self.setBcurvY:advance(0.0, 0.0, {}, {self.geo.bcurvY})
   else self.geo.bcurvY:clear(0.0) end
   if self.setPhiWall then self.setPhiWall:advance(0.0, 0.0, {}, {self.geo.phiWall})
   else self.geo.phiWall:clear(0.0) end
   -- sync ghost cells. these calls do not enforce periodicity because
   -- these fields initialized with syncPeriodicDirs = false
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.bcurvY:sync(false)
   self.geo.phiWall:sync(false)
end

function GkGeometry:write(tm)
   if self.evolve then
      if self.ioTrigger(tm) then
	 self.fieldIo:write(self.geo.bmag, string.format("bmag_%d.bp", self.ioFrame), tm)
	 self.fieldIo:write(self.geo.bcurvY, string.format("bcurvY_%d.bp", self.ioFrame), tm)
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.geo.bmag, string.format("bmag_%d.bp", self.ioFrame), tm)
	 self.fieldIo:write(self.geo.bcurvY, string.format("bcurvY_%d.bp", self.ioFrame), tm)
      end
      self.ioFrame = self.ioFrame+1
   end
end

function GkGeometry:rkStepperFields()
   return { self.geo, self.geo, self.geo, self.geo }
end

function GkGeometry:forwardEuler(tCurr, dt, momIn, geoIn, geoOut)
   if self.evolve then 
      self.setPhiWall:advance(tCurr, dt, {}, self.geo.phiWall)
   end 
   return true, GKYL_MAX_DOUBLE
end

function GkGeometry:applyBc(tCurr, dt, geoIn)
   if self.evolve then 
      self.geo.phiWall:sync(false)
   end
end

function GkGeometry:totalSolverTime()
   return self.setBmag.totalTime + self.setBcurvY.totalTime
end

function GkGeometry:totalBcTime() return 0.0 end
function GkGeometry:energyCalcTime() return 0.0 end


return {GkField = GkField, GkGeometry = GkGeometry}
