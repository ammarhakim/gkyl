-- GkField ---------------------------------------------------------------------
--
-- App support code: Gyrokinetic fields phi and apar, solved by
-- (perpendicular) Poisson and Ampere equations
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local DataStruct = require "DataStruct"
local LinearTrigger = require "LinearTrigger"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"
local FieldBase = require "App.Field.FieldBase"
local Species = require "App.Species"
local Time = require "Lib.Time"
local math = require("sci.math").generic
local diff = require("sci.diff")

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
   if tbl.phiBcBack then self.phiBcBack = tbl.phiBcBack end
   if tbl.phiBcFront then self.phiBcFront = tbl.phiBcFront end
   if tbl.aparBcLeft then self.aparBcLeft = tbl.aparBcLeft end
   if tbl.aparBcRight then self.aparBcRight = tbl.aparBcRight end
   if tbl.aparBcBottom then self.aparBcBottom = tbl.aparBcBottom end
   if tbl.aparBcTop then self.aparBcTop = tbl.aparBcTop end
   if appTbl.periodicDirs then self.periodicDirs = appTbl.periodicDirs end

   -- for storing integrated energies
   self.phi2 = DataStruct.DynVector { numComponents = 1 }
   self.apar2 = DataStruct.DynVector { numComponents = 1 }

   self.adiabatic = false
   self.discontinuousPhi = xsys.pickBool(tbl.discontinuousPhi, false)
   self.discontinuousApar = xsys.pickBool(tbl.discontinuousApar, true)

   -- for ndim=1 non-adiabatic only
   self.kperp2 = tbl.kperp2

   -- allow user to specify polarization weight. will be calculated automatically if not specified
   self.polarizationWeight = tbl.polarizationWeight

   if self.isElectromagnetic then
      self.mu0 = assert(tbl.mu0, "GkField: must specify mu0 for electromagnetic")
      self.dApardtInitFunc = tbl.dApardtInit
   end

   self.bcTime = 0.0 -- timer for BCs
end

-- methods for EM field object
function GkField:hasEB() return true, self.isElectromagnetic end
function GkField:setCfl() end
function GkField:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkField:setBasis(basis) self.basis = basis end
function GkField:setGrid(grid) self.grid = grid; self.ndim = self.grid:ndim() end

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
         self.potentials[i].dApardt = DataStruct.Field {
               onGrid = self.grid,
               numComponents = self.basis:numBasis(),
               ghost = {1, 1}
         }
         self.potentials[i].apar:clear(0.0)
         self.potentials[i].dApardt:clear(0.0)
      end
   end

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
   self.massWeight = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
   }
end

-- solve for initial fields self-consistently 
-- from initial distribution function
function GkField:initField(species)
   -- solve for initial phi
   self:forwardEuler(0, 1.0, species, 1, 1)

   if self.isElectromagnetic then
      -- solve for initial Apar
      self.currentDens:clear(0.0)
      for nm, s in pairs(species) do
         self.currentDens:accumulate(s:getCharge(), s:getMomDensity())
      end
      self.aparSlvr:advance(tCurr, dt, {self.currentDens}, {self.potentials[1].apar})

      -- clear dApar/dt ... will be solved for before being used
      self.potentials[1].dApardt:clear(0.0)
   end

   -- apply BCs 
   self:applyBc(0, 0, self.potentials[1])
end

function GkField:rkStepperFields()
   return self.potentials
end

-- for RK timestepping for non-elliptic fields (e.g. only apar)
function GkField:copyRk(outIdx, aIdx)
   if self.isElectromagnetic and self:rkStepperFields()[aIdx] then 
      self:rkStepperFields()[outIdx].apar:copy(self:rkStepperFields()[aIdx].apar)
   end
end
-- for RK timestepping for non-elliptic fields (e.g. only apar)
function GkField:combineRk(outIdx, a, aIdx, ...)
   if self.isElectromagnetic and self:rkStepperFields()[aIdx] then
      local args = {...} -- package up rest of args as table
      local nFlds = #args/2
      self:rkStepperFields()[outIdx].apar:combine(a, self:rkStepperFields()[aIdx].apar)
      for i = 1, nFlds do -- accumulate rest of the fields
         self:rkStepperFields()[outIdx].apar:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]].apar)
      end	 
   end
end

function GkField:createSolver(species)
   -- leaving this here for a possible future implementation...
   --local laplacianWeight = {}
   --local modifierConstant = {}
   --if self.adiabatic then
   --  local fac = self.adiabQneutFac
   --  for d = 1, self.ndim do 
   --    laplacianWeight[d] = 0.0
   --    modifierConstant[d] = fac
   --  end
   --else
   --  if ndim==1 then  -- z
   --    laplacianWeight = {0.0}
   --    modifierConstant = {self.kperp2}
   --  elseif ndim==2 then  -- x,y
   --    laplacianWeight = {1.0, 1.0}
   --    modifierConstant = {0.0, 0.0}
   --  else -- x,y,z
   --    laplacianWeight = {1.0, 1.0, 0.0}
   --    modifierConstant = {0.0, 0.0, 1.0}
   --  end
   --end

   -- get adiabatic species info and calculate species-dependent 
   -- weight on polarization term == sum_s m_s n_s / B^2
   if not self.polarizationWeight then 
      self.polarizationWeight = 0.0
      for nm, s in pairs(species) do
         if Species.AdiabaticSpecies.is(s) then
            self.adiabatic = true
            self.adiabSpec = s
         elseif Species.GkSpecies.is(s) then
            self.polarizationWeight = self.polarizationWeight + s:polarizationWeight()
         end
      end
   end
   assert((self.adiabatic and self.isElectromagnetic) == false, "GkField: cannot use adiabatic response for electromagnetic case")

   if self.adiabatic then
      -- if adiabatic, we don't solve the Poisson equation, but
      -- we do need to project the discontinuous charge densities
      -- onto a continuous basis to set phi. 
      self.phiSlvr = Updater.FemPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcLeft = self.phiBcLeft,
        bcRight = self.phiBcRight,
        bcBottom = self.phiBcBottom,
        bcTop = self.phiBcTop,
        bcBack = self.phiBcBack,
        bcFront = self.phiBcFront,
        periodicDirs = self.periodicDirs,
        laplacianWeight = 0.0, 
        modifierConstant = self.adiabSpec:qneutFac(), -- = n0*q^2/T
        zContinuous = not self.discontinuousPhi,
      }
   else
      local ndim = self.ndim
      local laplacianWeight, modifierConstant
      assert(self.polarizationWeight, "GkField: must specify polarizationWeight = ni*mi/B^2 for non-adiabatic field")
      if ndim==1 then  -- z
        assert(self.kperp2, "GkField: must specify kperp2 for non-adiabatic field with ndim=1")
        laplacianWeight = 0.0 -- {0.0}
        modifierConstant = self.kperp2*self.polarizationWeight -- {self.kperp2}
      elseif ndim==2 then  -- x,y
        laplacianWeight = self.polarizationWeight -- {1.0, 1.0}
        modifierConstant = 0.0 --  {0.0, 0.0}
      else -- x,y,z
        laplacianWeight = self.polarizationWeight -- {1.0, 1.0, 0.0}
        modifierConstant = 0.0 -- {0.0, 0.0, 1.0}
      end
      self.phiSlvr = Updater.FemPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcLeft = self.phiBcLeft,
        bcRight = self.phiBcRight,
        bcBottom = self.phiBcBottom,
        bcTop = self.phiBcTop,
        periodicDirs = self.periodicDirs,
        laplacianWeight = laplacianWeight,
        modifierConstant = modifierConstant,
        zContinuous = not self.discontinuousPhi,
      }
      if self.isElectromagnetic then 
        -- NOTE: aparSlvr only used to solve for initial Apar
        -- at all other times Apar is timestepped using dApar/dt
        if ndim==1 then
           laplacianWeight = 0.0
           modifierConstant = 1.0/self.mu0
        else
           laplacianWeight = 1.0/self.mu0
           modifierConstant = 0.0
        end
        self.aparSlvr = Updater.FemPoisson {
          onGrid = self.grid,
          basis = self.basis,
          bcLeft = self.aparBcLeft,
          bcRight = self.aparBcRight,
          bcBottom = self.aparBcBottom,
          bcTop = self.aparBcTop,
          periodicDirs = self.periodicDirs,
          laplacianWeight = laplacianWeight,
          modifierConstant = modifierConstant,
          zContinuous = not self.discontinuousApar,
        }

        if ndim==1 then
           laplacianWeight = 0.0
           modifierConstant = 1.0
        else
           laplacianWeight = 1.0/self.mu0
           modifierConstant = 1.0
        end
        self.dApardtSlvr = Updater.FemPoisson {
          onGrid = self.grid,
          basis = self.basis,
          bcLeft = self.aparBcLeft,
          bcRight = self.aparBcRight,
          bcBottom = self.aparBcBottom,
          bcTop = self.aparBcTop,
          periodicDirs = self.periodicDirs,
          laplacianWeight = laplacianWeight,
          modifierConstant = modifierConstant,
          zContinuous = not self.discontinuousApar,
        }

        -- set up constant dummy field
        self.unitWeight = DataStruct.Field {
             onGrid = self.grid,
             numComponents = self.basis:numBasis(),
             ghost = {1, 1},
        }
        local initUnit = Updater.ProjectOnBasis {
           onGrid = self.grid,
           basis = self.basis,
           evaluate = function (t,xn)
                         return 1.0
                      end
        }
        initUnit:advance(0.,0.,{},{self.unitWeight})
      end
   end

   -- need to set this flag so that field calculated self-consistently at end of full RK timestep
   self.isElliptic = true

   if self.isElectromagnetic then self.step2 = true end
end

function GkField:createDiagnostics()
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.potentials[1].phi:elemType(),
      method = self.ioMethod,
   }

   self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      quantity = "V2"
   }
end

function GkField:write(tm)
   if self.evolve then
      -- compute integrated quantities over domain
      self.energyCalc:advance(tm, 0.0, { self.potentials[1].phi }, { self.phi2 })
      if self.isElectromagnetic then 
        self.energyCalc:advance(tm, 0.0, { self.potentials[1].apar }, { self.apar2 })
      end
      
      if self.ioTrigger(tm) then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
	   self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
	 self.phi2:write(string.format("phi2_%d.bp", self.ioFrame), tm, self.ioFrame)
	 if self.isElectromagnetic then
	    self.apar2:write(string.format("apar2_%d.bp", self.ioFrame), tm, self.ioFrame)
	 end
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
	   self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
      end
      self.ioFrame = self.ioFrame+1
   end
end

function GkField:writeRestart(tm)
end

-- not needed for GK
function GkField:accumulateCurrent(dt, current, em)
end

function GkField:forwardEuler(tCurr, dt, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potNew = self:rkStepperFields()[outIdx]

   if self.evolve or tCurr == 0.0 then
      local mys, mys2 = true, true
      local mydt, mydt2 = GKYL_MAX_DOUBLE, GKYL_MAX_DOUBLE
      self.chargeDens:clear(0.0)
      for nm, s in pairs(species) do
         self.chargeDens:accumulate(s:getCharge(), s:getNumDensity())
      end
      -- phi solve (elliptic, so update potCurr.phi)
      local mys, mydt = self.phiSlvr:advance(tCurr, dt, {self.chargeDens}, {potCurr.phi})

      -- apply BCs
      local tmStart = Time.clock()
      potCurr.phi:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)
 
      return mys and mys2, math.min(mydt,mydt2)
   else
      -- just copy stuff over
      potNew.phi:copy(potCurr.phi)
      if self.isElectromagnetic then 
         potNew.apar:copy(potCurr.apar) 
         potNew.dApardt:copy(potCurr.dApardt) 
      end
      return true, GKYL_MAX_DOUBLE
   end
end

function GkField:forwardEulerStep2(tCurr, dt, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potNew = self:rkStepperFields()[outIdx]

   if self.evolve then
      local mys, mys2 = true, true
      local mydt, mydt2 = GKYL_MAX_DOUBLE, GKYL_MAX_DOUBLE
      self.currentDens:clear(0.0)
      if self.ndim==1 then 
         self.massWeight:combine(self.kperp2/self.mu0, self.unitWeight) 
      else 
         self.massWeight:clear(0.0)
      end
      for nm, s in pairs(species) do
        if s:isEvolving() then 
           self.massWeight:accumulate(s:getCharge()*s:getCharge()/s:getMass(), s:getNumDensity())
           -- taking momDensity at outIdx gives momentum moment of df/dt
           self.currentDens:accumulate(s:getCharge(), s:getMomDensity(outIdx))
        end
      end
      -- dApar/dt solve
      local dApardt = potCurr.dApardt
      mys2, mydt2 = self.dApardtSlvr:advance(tCurr, dt, {self.currentDens, nil, self.massWeight}, {dApardt}) 
      
      -- decrease effective polynomial order in z of dApar/dt by setting the highest order z coefficients to 0
      -- this ensures that dApar/dt is in the same space as dPhi/dz
      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
         local polyOrder = self.basis:polyOrder()
         local localRange = dApardt:localRange()
         local indexer = dApardt:genIndexer()
         local ptr = dApardt:get(1)

         -- loop over all cells
         for idx in localRange:colMajorIter() do
            self.grid:setIndex(idx)
            
            dApardt:fill(indexer(idx), ptr)
            if self.ndim == 1 then
               ptr:data()[polyOrder] = 0.0
            else -- ndim == 3
               if polyOrder == 1 then
                  ptr:data()[3] = 0.0
                  ptr:data()[5] = 0.0
                  ptr:data()[6] = 0.0
                  ptr:data()[7] = 0.0
               elseif polyOrder == 2 then
                  ptr:data()[9] = 0.0
                  ptr:data()[13] = 0.0
                  ptr:data()[14] = 0.0
                  ptr:data()[15] = 0.0
                  ptr:data()[16] = 0.0
                  ptr:data()[17] = 0.0
                  ptr:data()[18] = 0.0
                  ptr:data()[19] = 0.0
               end
            end
         end
      end

      -- advance Apar
      potNew.apar:combine(1.0, potCurr.apar, dt, dApardt)

      -- apply BCs
      local tmStart = Time.clock()
      dApardt:sync(true)
      potNew.apar:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      return mys and mys2, math.min(mydt,mydt2)
   else
      return true, GKYL_MAX_DOUBLE
   end
end

-- NOTE: global boundary conditions handled by solver. this just updates interproc ghosts.
function GkField:applyBc(tCurr, dt, potIn)
   local tmStart = Time.clock()
   potIn.phi:sync(true)
   if self.isElectromagnetic then 
     potIn.apar:sync(true) 
     potIn.dApardt:sync(true) 
   end
   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end
   
function GkField:totalSolverTime()
   if self.phiSlvr then
     local time = self.phiSlvr.totalTime
     if self.isElectromagnetic and self.aparSlvr then time = time + self.aparSlvr.totalTime end
     if self.isElectromagnetic and self.dApardtSlvr then time = time + self.dApardtSlvr.totalTime end
     return time
   end
   return 0.0
end

function GkField:totalBcTime()
   return self.bcTime
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

   -- scalar magnetic field strength (e.g. to calculate gyro-radii with)
   -- calculate from max(bmag)? 
   self.B0 = assert(B0 or tbl.B0, "Must specify B0 as global App variable or in GkGeometry as 'B0'")
   
   -- get function to initialize background magnetic field
   self.bmagFunc = tbl.bmag
   assert(self.bmagFunc and type(self.bmagFunc)=="function", "GkGeometry: must specify background magnetic field function with 'bmag'")

   -- wall potential for sheath BCs
   self.phiWallFunc = tbl.phiWall
   if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "GkGeometry: phiWall must be a function (t, xn)") end
end

function GkGeometry:hasEB() end
function GkGeometry:setCfl() end
function GkGeometry:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkGeometry:setBasis(basis) self.basis = basis end
function GkGeometry:setGrid(grid) self.grid = grid; self.ndim = self.grid:ndim() end

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

   -- functions for magnetic drifts ~ 1/B*curl(bhat) 
   self.geo.bdriftX = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   self.geo.bdriftY = DataStruct.Field {
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
   -- determine which variables bmag depends on by checking if setting a variable to nan results in nan
   local ones = {}
   for dir = 1, self.ndim do
      ones[dir] = 1
   end
   self.bmagVars = {}
   for dir = 1, self.ndim do
      ones[dir] = 0/0 -- set this var to nan 
      -- test if result is nan.. nan is the only value that doesn't equal itself
      if self.bmagFunc(0, ones) ~= self.bmagFunc(0, ones) then 
        -- if result is nan, bmag must depend on this var
        table.insert(self.bmagVars, dir) 
      end
      ones[dir] = 1 -- reset so we can check other vars
   end
   if self.bmagVars[1] == nil then self.bmagVars[1] = 0 end

   -- calculate 1/B function
   self.bmagInvFunc = function (t, xn)
      return 1/self.bmagFunc(t,xn)
   end
   -- calculate magnetic drift functions
   if self.ndim > 1 then
      local function bgrad(xn)
         local function bmagUnpack(...)
            local xn = {...}
            return self.bmagFunc(0, xn)
         end
         local deriv = diff.derivativef(bmagUnpack, #xn)
         xntable = {}
         for i = 1, #xn do
           xntable[i] = xn[i]
         end
         local f, dx, dy, dz = deriv(unpack(xntable))
         return dx, dy, dz
      end
      self.bdriftXFunc = function (t, xn)
         local bgradX, bgradY, bgradZ = bgrad(xn)
         return -bgradY/self.bmagFunc(t,xn)^2
      end
      self.bdriftYFunc = function (t, xn)
         local bgradX, bgradY, bgradZ = bgrad(xn)
         return bgradX/self.bmagFunc(t,xn)^2
      end
   end

   -- projection updaters
   self.setBmag = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.bmagFunc,
      projectOnGhosts = true,
   }
   self.setBmagInv = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.bmagInvFunc
   }
   if self.bdriftXFunc then 
      self.setBdriftX = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         projectOnGhosts = true,
         evaluate = self.bdriftXFunc
      }
   end
   if self.bdriftYFunc then
      self.setBdriftY = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         projectOnGhosts = true,
         evaluate = self.bdriftYFunc
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
   if self.setBdriftX then self.setBdriftX:advance(0.0, 0.0, {}, {self.geo.bdriftX})
   else self.geo.bdriftX:clear(0.0) end
   if self.setBdriftY then self.setBdriftY:advance(0.0, 0.0, {}, {self.geo.bdriftY})
   else self.geo.bdriftY:clear(0.0) end
   if self.setPhiWall then self.setPhiWall:advance(0.0, 0.0, {}, {self.geo.phiWall})
   else self.geo.phiWall:clear(0.0) end
   -- sync ghost cells. these calls do not enforce periodicity because
   -- these fields initialized with syncPeriodicDirs = false
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.bdriftX:sync(false)
   self.geo.bdriftY:sync(false)
   self.geo.phiWall:sync(false)

   -- apply BCs
   self:applyBc(0, 0, self.geo)
end

function GkGeometry:write(tm)
   -- not evolving geometry, so only write geometry at beginning
   if self.ioFrame == 0 then
      self.fieldIo:write(self.geo.bmag, string.format("bmag_%d.bp", self.ioFrame), tm, self.ioFrame)
      if self.setBdriftX then
	 self.fieldIo:write(self.geo.bdriftX, string.format("bdriftX_%d.bp", self.ioFrame), tm, self.ioFrame)
      end
      if self.setBdriftY then
	 self.fieldIo:write(self.geo.bdriftY, string.format("bdriftY_%d.bp", self.ioFrame), tm, self.ioFrame)
      end
   end
   self.ioFrame = self.ioFrame+1
end

function GkGeometry:writeRestart(tm)
end

function GkGeometry:rkStepperFields()
   return { self.geo, self.geo, self.geo, self.geo }
end

function GkGeometry:forwardEuler(tCurr, dt)
   if self.evolve then 
      self.setPhiWall:advance(tCurr, dt, {}, self.geo.phiWall)
   end 
   return true, GKYL_MAX_DOUBLE
end

function GkGeometry:applyBc(tCurr, dt, geoIn)
   if self.evolve then 
      geoIn.phiWall:sync(false)
   end
end

function GkGeometry:totalSolverTime()
   if self.evolve then return self.setPhiWall.totalTime
   else return 0 end
end

function GkGeometry:totalBcTime() return 0.0 end
function GkGeometry:energyCalcTime() return 0.0 end


return {GkField = GkField, GkGeometry = GkGeometry}
