-- GkField ---------------------------------------------------------------------
--
-- App support code: Gyrokinetic fields phi and apar, solved by
-- (perpendicular) Poisson and Ampere equations
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Constants = require "Lib.Constants"
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
local Constants = require "Lib.Constants"

local GkField = Proto(FieldBase.FieldBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function GkField:init(tbl)
   GkField.super.init(self, tbl)
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

   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

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
   self.esEnergy = DataStruct.DynVector { numComponents = 1 }
   self.emEnergy = DataStruct.DynVector { numComponents = 1 }

   self.adiabatic = false
   self.discontinuousPhi = xsys.pickBool(tbl.discontinuousPhi, false)
   self.discontinuousApar = xsys.pickBool(tbl.discontinuousApar, true)

   -- for ndim=1 only
   self.kperp2 = tbl.kperp2

   -- allow user to specify polarization weight. will be calculated automatically if not specified
   self.polarizationWeight = tbl.polarizationWeight

   -- determine whether to use linearized polarization term in poisson equation, which uses background density in polarization weight
   -- if not, uses full time-dependent density in polarization weight 
   self.linearizedPolarization = xsys.pickBool(tbl.linearizedPolarization, true)
   self.uniformPolarization = xsys.pickBool(tbl.uniformPolarization, true)

   if self.isElectromagnetic then
      self.mu0 = tbl.mu0 or Constants.MU0
   end

   self.initPhiFunc = tbl.initPhiFunc
   if self.initPhiFunc and self.evolve then 
      print("GkField: warning... specifying initPhiFunc will make initial phi inconsistent with f") 
   end

   self.bcTime = 0.0 -- timer for BCs

   self._first = true
   self._firstStep = true
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

   self.dApardtProv = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
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
                 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- set up some other fields
   self.weight = DataStruct.Field {
        onGrid = self.grid,
        numComponents = self.basis:numBasis(),
        ghost = {1, 1},
   }
   self.laplacianWeight = DataStruct.Field {
        onGrid = self.grid,
        numComponents = self.basis:numBasis(),
        ghost = {1, 1},
   }
   self.modifierWeight = DataStruct.Field {
        onGrid = self.grid,
        numComponents = self.basis:numBasis(),
        ghost = {1, 1},
   }
end

-- solve for initial fields self-consistently 
-- from initial distribution function
function GkField:initField(species)
   if self.initPhiFunc then
      local project = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self.initPhiFunc,
         projectOnGhosts = true
      }
      project:advance(0.0, {}, {self.potentials[1].phi})
   else
      -- solve for initial phi
      self:advance(0.0, species, 1, 1)
   end

   if self.isElectromagnetic then
      -- solve for initial Apar
      self.currentDens:clear(0.0)
      for nm, s in pairs(species) do
         self.currentDens:accumulate(s:getCharge(), s:getMomDensity())
      end
      self.aparSlvr:advance(0.0, {self.currentDens}, {self.potentials[1].apar})

      -- clear dApar/dt ... will be solved for before being used
      self.potentials[1].dApardt:clear(0.0)
   end

   -- apply BCs 
   self:applyBc(0, self.potentials[1])
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

function GkField:suggestDt()
   return GKYL_MAX_DOUBLE
end
function GkField:clearCFL()
end

function GkField:createSolver(species, funcField)
   -- get adiabatic species info
   for nm, s in pairs(species) do
      if Species.AdiabaticSpecies.is(s) then
         self.adiabatic = true
         self.adiabSpec = s
      end
   end
   assert((self.adiabatic and self.isElectromagnetic) == false, "GkField: cannot use adiabatic response for electromagnetic case")

   -- set up FEM solver for Poisson equation to solve for phi
   local gxx, gxy, gyy, jacobGeo
   if funcField.geo then 
     gxx = funcField.geo.gxx
     gxy = funcField.geo.gxy
     gyy = funcField.geo.gyy
     jacobGeo = funcField.geo.jacobGeo
   end
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
     zContinuous = not self.discontinuousPhi,
     gxx = gxx,
     gxy = gxy,
     gyy = gyy,
   }
   -- when using a linearizedPolarization term in Poisson equation,
   -- the weights on the terms are constant scalars
   if self.linearizedPolarization then
      -- if not provided, calculate species-dependent weight on polarization term == sum_s m_s n_s / B^2
      if not self.polarizationWeight then 
         self.polarizationWeight = 0.0
         for nm, s in pairs(species) do
            if Species.GkSpecies.is(s) then
               self.polarizationWeight = self.polarizationWeight + s:getPolarizationWeight()
            end
         end
      end
      -- if not adiabatic, and polarization weight still not set, assume it is 1
      if self.polarizationWeight == 0.0 and not self.adiabatic then self.polarizationWeight = 1.0 end

      -- set up scalar multipliers on laplacian and modifier terms in Poisson equation
      local laplacianConstant, modifierConstant
 
      if self.ndim==1 then
         assert(self.kperp2, "GkField: must specify kperp2 for ndim=1")
         laplacianConstant = 0.0 
         modifierConstant = self.kperp2*self.polarizationWeight 
      else 
         laplacianConstant = -self.polarizationWeight 
         modifierConstant = 0.0 
      end

      if self.adiabatic then 
         modifierConstant = modifierConstant + self.adiabSpec:getQneutFac() 
      end

      self.laplacianWeight:combine(laplacianConstant, jacobGeo)
      self.modifierWeight:combine(modifierConstant, jacobGeo)

      if laplacianConstant ~= 0 then self.phiSlvr:setLaplacianWeight(self.laplacianWeight) end
      if modifierConstant ~= 0 then self.phiSlvr:setModifierWeight(self.modifierWeight) end
   end
   -- when not using linearizedPolarization, weights are set each step in advance method

   if self.isElectromagnetic then 
     local ndim = self.ndim
     local laplacianConstant, modifierConstant
     -- NOTE: aparSlvr only used to solve for initial Apar
     -- at all other times Apar is timestepped using dApar/dt
     self.aparSlvr = Updater.FemPoisson {
       onGrid = self.grid,
       basis = self.basis,
       bcLeft = self.aparBcLeft,
       bcRight = self.aparBcRight,
       bcBottom = self.aparBcBottom,
       bcTop = self.aparBcTop,
       periodicDirs = self.periodicDirs,
       zContinuous = not self.discontinuousApar,
       gxx = gxx,
       gxy = gxy,
       gyy = gyy,
     }
     if ndim==1 then
        laplacianConstant = 0.0
        modifierConstant = self.kperp2/self.mu0
     else
        laplacianConstant = -1.0/self.mu0
        modifierConstant = 0.0
     end
     self.laplacianWeight:combine(laplacianConstant, self.unitWeight)
     self.modifierWeight:combine(modifierConstant, self.unitWeight)
     if laplacianConstant ~= 0 then self.aparSlvr:setLaplacianWeight(self.laplacianWeight) end
     if modifierConstant ~= 0 then self.aparSlvr:setModifierWeight(self.modifierWeight) end

     self.dApardtSlvr = Updater.FemPoisson {
       onGrid = self.grid,
       basis = self.basis,
       bcLeft = self.aparBcLeft,
       bcRight = self.aparBcRight,
       bcBottom = self.aparBcBottom,
       bcTop = self.aparBcTop,
       periodicDirs = self.periodicDirs,
       zContinuous = not self.discontinuousApar,
       gxx = gxx,
       gxy = gxy,
       gyy = gyy,
     }
     if ndim==1 then
        laplacianConstant = 0.0
        modifierConstant = 1.0
     else
        laplacianConstant = -1.0/self.mu0
        modifierConstant = 1.0
     end
     self.laplacianWeight:combine(laplacianConstant, self.unitWeight)
     self.modifierWeight:combine(modifierConstant, self.unitWeight)
     if laplacianConstant ~= 0 then self.dApardtSlvr:setLaplacianWeight(self.laplacianWeight) end
     if modifierConstant ~= 0 then self.dApardtSlvr:setModifierWeight(self.modifierWeight) end
   end

   -- need to set this flag so that field calculated self-consistently at end of full RK timestep
   self.isElliptic = true

   if self.isElectromagnetic and self.basis:polyOrder() == 1 then 
      self.nstep = 3 
   elseif self.isElectromagnetic then
      self.nstep = 2
   end
end

function GkField:createDiagnostics()
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.potentials[1].phi:elemType(),
      method = self.ioMethod,
      writeGhost = self.writeGhost
   }

   -- updaters for computing integrated quantities
   self.int2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      quantity = "V2"
   }
   if self.ndim == 1 then
      self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,
         basis = self.basis,
         quantity = "V2"
      }
   elseif self.basis:polyOrder()==1 then -- GradPerpV2 only implemented for p=1 currently
      self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,
         basis = self.basis,
         quantity = "GradPerpV2"
      }
   end
end

function GkField:write(tm, force)
   if self.evolve then
      -- compute integrated quantities over domain
      self.int2Calc:advance(tm, { self.potentials[1].phi }, { self.phi2 })
      if self.isElectromagnetic then 
        self.int2Calc:advance(tm, { self.potentials[1].apar }, { self.apar2 })
      end
      local esEnergyFac = .5*self.polarizationWeight
      if self.ndim == 1 then esEnergyFac = esEnergyFac*self.kperp2 end
      if self.energyCalc then self.energyCalc:advance(tm, { self.potentials[1].phi, esEnergyFac }, { self.esEnergy }) end
      if self.isElectromagnetic then 
        local emEnergyFac = .5/self.mu0
        if self.ndim == 1 then emEnergyFac = emEnergyFac*self.kperp2 end
        if self.energyCalc then self.energyCalc:advance(tm, { self.potentials[1].apar, emEnergyFac}, { self.emEnergy }) end
      end
      
      if self.ioTrigger(tm) or force then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
	   self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
	 self.phi2:write(string.format("phi2_%d.bp", self.ioFrame), tm, self.ioFrame)
	 self.esEnergy:write(string.format("esEnergy_%d.bp", self.ioFrame), tm, self.ioFrame)
	 if self.isElectromagnetic then
	    self.apar2:write(string.format("apar2_%d.bp", self.ioFrame), tm, self.ioFrame)
	    self.emEnergy:write(string.format("emEnergy_%d.bp", self.ioFrame), tm, self.ioFrame)
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
   -- (the final "false" prevents writing of ghost cells)
   self.fieldIo:write(self.potentials[1].phi, "phi_restart.bp", tm, self.ioFrame, false)
   self.fieldIo:write(self.phiSlvr:getLaplacianWeight(), "laplacianWeight_restart.bp", tm, self.ioFrame, false)
   self.fieldIo:write(self.phiSlvr:getModifierWeight(), "modifierWeight_restart.bp", tm, self.ioFrame, false)
   if self.isElectromagnetic then
     self.fieldIo:write(self.potentials[1].apar, "apar_restart.bp", tm, self.ioFrame, false)
     self.fieldIo:write(self.potentials[1].dApardt, "dApardt_restart.bp", tm, self.ioFrame, false)
     self.fieldIo:write(self.dApardtSlvr:getModifierWeight(), "modifierWeightEM_restart.bp", tm, self.ioFrame, false)
   end

   -- (the final "false" prevents flushing of data after write)
   self.phi2:write("phi2_restart.bp", tm, self.ioFrame, false)

end

function GkField:readRestart()
   -- this read of restart file of phi is only to get frame
   -- numbering correct. The forward Euler recomputes the potential
   -- before updating the hyperbolic part
   local tm, fr = self.fieldIo:read(self.potentials[1].phi, "phi_restart.bp")
   self.phi2:read("phi2_restart.bp", tm)

   self.fieldIo:read(self.laplacianWeight, "laplacianWeight_restart.bp")
   self.fieldIo:read(self.modifierWeight, "modifierWeight_restart.bp")
   self.phiSlvr:setLaplacianWeight(self.laplacianWeight)
   if self.adiabatic or self.ndim == 1 then self.phiSlvr:setModifierWeight(self.modifierWeight) end

   if self.isElectromagnetic then
      self.fieldIo:read(self.potentials[1].apar, "apar_restart.bp")
      self.fieldIo:read(self.potentials[1].dApardt, "dApardt_restart.bp")
      self.fieldIo:read(self.modifierWeight, "modifierWeightEM_restart.bp")
      self.dApardtSlvr:setModifierWeight(self.modifierWeight)
   end

   self:applyBc(0, self.potentials[1])

   self.ioFrame = fr 
   -- iterate triggers
   self.ioTrigger(tm)
end

-- not needed for GK
function GkField:accumulateCurrent(dt, current, em)
end

-- solve for electrostatic potential phi
function GkField:advance(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs = self:rkStepperFields()[outIdx]
   
   if self.evolve or (self._first and not self.initPhiFunc) then
      self.chargeDens:clear(0.0)
      for nm, s in pairs(species) do
         self.chargeDens:accumulate(s:getCharge(), s:getNumDensity())
      end
      -- if not using linearized polarization term, set up laplacian weight
      if not self.linearizedPolarization or (self._first and not self.uniformPolarization) then
         self.weight:clear(0.0)
         for nm, s in pairs(species) do
            if Species.GkSpecies.is(s) then
               self.weight:accumulate(1.0, s:getPolarizationWeight(false))
            end
         end
         if self.ndim == 1 then
            self.modifierWeight:combine(self.kperp2, self.weight)
         else
            self.modifierWeight:clear(0.0)
            self.laplacianWeight:combine(-1.0, self.weight)
         end

         if self.adiabatic then
            self.modifierWeight:accumulate(1.0, self.adiabSpec:getQneutFac(false))
         end

         if self.ndim > 1 then
            self.phiSlvr:setLaplacianWeight(self.laplacianWeight)
         end
         if self.adiabatic or self.ndim == 1 then self.phiSlvr:setModifierWeight(self.modifierWeight) end
      end
      -- phi solve (elliptic, so update potCurr.phi)
      self.phiSlvr:advance(tCurr, {self.chargeDens}, {potCurr.phi})

      -- apply BCs
      local tmStart = Time.clock()
      potCurr.phi:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)
 
      self._first = false
   else
      -- just copy stuff over
      if self.isElectromagnetic then 
         potRhs.apar:copy(potCurr.apar) 
      end
   end
end

-- solve for dApardt in p>=2, or solve for a provisional dApardtProv in p=1
function GkField:advanceStep2(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs = self:rkStepperFields()[outIdx]

   local polyOrder = self.basis:polyOrder()

   if self.evolve then

      self.currentDens:clear(0.0)
      if self.ndim==1 then 
         self.modifierWeight:combine(self.kperp2/self.mu0, self.unitWeight) 
      else 
         self.modifierWeight:clear(0.0)
      end
      for nm, s in pairs(species) do
        if s:isEvolving() then 
           self.modifierWeight:accumulate(s:getCharge()*s:getCharge()/s:getMass(), s:getNumDensity())
           -- taking momDensity at outIdx gives momentum moment of df/dt
           self.currentDens:accumulate(s:getCharge(), s:getMomDensity(outIdx))
        end
      end
      self.dApardtSlvr:setModifierWeight(self.modifierWeight)
      -- dApar/dt solve
      local dApardt = potCurr.dApardt
      self.dApardtSlvr:advance(tCurr, {self.currentDens}, {dApardt}) 
      
      -- decrease effective polynomial order in z of dApar/dt by setting the highest order z coefficients to 0
      -- this ensures that dApar/dt is in the same space as dPhi/dz
      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
         local localRange = dApardt:localRange()
         local indexer = dApardt:genIndexer()
         local ptr = dApardt:get(1)

         -- loop over all cells
         for idx in localRange:rowMajorIter() do
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

      -- apply BCs
      local tmStart = Time.clock()
      dApardt:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      if polyOrder > 1 then 
         -- Apar RHS is just dApar/dt
         potRhs.apar:copy(dApardt)
      else 
         -- save dApardt as dApardtProv, so that it can be used in upwinding 
         -- in vpar surface terms in p=1 Ohm's law and GK update
         self.dApardtProv:copy(dApardt)
      end
   end
end

-- note: step 3 is for p=1 only: solve for dApardt
function GkField:advanceStep3(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs = self:rkStepperFields()[outIdx]

   local polyOrder = self.basis:polyOrder()

   if self.evolve then
      self.currentDens:clear(0.0)
      if self.ndim==1 then 
         self.modifierWeight:combine(self.kperp2/self.mu0, self.unitWeight) 
      else 
         self.modifierWeight:clear(0.0)
      end
      for nm, s in pairs(species) do
        if s:isEvolving() then 
           self.modifierWeight:accumulate(s:getCharge()*s:getCharge()/s:getMass(), s:getEmModifier())
           -- taking momDensity at outIdx gives momentum moment of df/dt
           self.currentDens:accumulate(s:getCharge(), s:getMomProjDensity(outIdx))
        end
      end
      self.dApardtSlvr:setModifierWeight(self.modifierWeight)
      -- dApar/dt solve
      local dApardt = potCurr.dApardt
      self.dApardtSlvr:advance(tCurr, {self.currentDens}, {dApardt}) 
      
      -- decrease effective polynomial order in z of dApar/dt by setting the highest order z coefficients to 0
      -- this ensures that dApar/dt is in the same space as dPhi/dz
      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
         local localRange = dApardt:localRange()
         local indexer = dApardt:genIndexer()
         local ptr = dApardt:get(1)

         -- loop over all cells
         for idx in localRange:rowMajorIter() do
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

      -- apply BCs
      local tmStart = Time.clock()
      dApardt:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      -- Apar RHS is just dApar/dt
      potRhs.apar:copy(dApardt)
   end
end

-- NOTE: global boundary conditions handled by solver. this just updates interproc ghosts.
function GkField:applyBc(tCurr, potIn)
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
   local t = self.int2Calc.totalTime
   if self.energyCalc then t = t + self.energyCalc.totalTime end
   return t
end

-- GkGeometry ---------------------------------------------------------------------
--
-- A field object with fields specifying the magnetic geometry for GK
--------------------------------------------------------------------------------

local GkGeometry = Proto(FieldBase.FuncFieldBase)

-- methods for no field object
function GkGeometry:init(tbl)
   GkGeometry.super.init(self, tbl)
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
   self.bmagFunc = tbl.bmag
   --assert(self.bmagFunc and type(self.bmagFunc)=="function", "GkGeometry: must specify background magnetic field function with 'bmag'")

   -- for s-alpha geometry
   self.salpha = tbl.salpha

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

   -- bmagInv = 1/B
   self.geo.bmagInv = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- bhat.grad z = 1/sqrt(g_zz)
   self.geo.gradpar = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- jacobian of coordinate transformation
   self.geo.jacobGeo = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- inverse of jacobian of coordinate transformation
   self.geo.jacobGeoInv = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- functions for magnetic drifts 
   -- geoX = g_xz / ( J B sqrt(g_zz) )
   self.geo.geoX = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- geoY = g_yz / ( J B sqrt(g_zz) )
   self.geo.geoY = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- geoZ = sqrt(g_zz) / ( J B )
   self.geo.geoZ = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
 
   -- functions for laplacian
   -- g^xx = |grad x|**2
   self.geo.gxx = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- g^xy = grad x . grad y
   self.geo.gxy = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- g^yy = |grad y|**2
   self.geo.gyy = DataStruct.Field {
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

   -- get needed metric functions from grid or defaults
   if self.ndim == 1 then
      self.g_zzFunc = function (t, xn)
         local g = {}
         self.grid:calcMetric(xn, g)
         return g[1]
      end
      self.g_xzFunc = function (t, xn)
         return 0
      end
      self.g_yzFunc = function (t, xn)
         return 0
      end
      self.gxx_Func = function (t, xn)
         return 1
      end
      self.gxy_Func = function (t, xn)
         return 0
      end
      self.gyy_Func = function (t, xn)
         return 1
      end
      self.jacobGeoFunc = function (t, xn)
         return self.grid:calcJacobian(xn)
      end
   elseif self.ndim == 2 then
      self.g_zzFunc = function (t, xn)
         return 1
      end
      self.g_xzFunc = function (t, xn)
         return 0
      end
      self.g_yzFunc = function (t, xn)
         return 0
      end
      self.gxx_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[1]
      end
      self.gxy_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[2]
      end
      self.gyy_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[3]
      end
      self.jacobGeoFunc = function (t, xn)
         return self.grid:calcJacobian(xn)
      end
   elseif self.ndim == 3 then
      self.g_zzFunc = function (t, xn)
         local g = {}
         self.grid:calcMetric(xn, g)
         return g[6]
      end
      self.g_xzFunc = function (t, xn)
         local g = {}
         self.grid:calcMetric(xn, g)
         return g[3]
      end
      self.g_yzFunc = function (t, xn)
         local g = {}
         self.grid:calcMetric(xn, g)
         return g[5]
      end
      self.gxx_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[1]
      end
      self.gxy_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[2]
      end
      self.gyy_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[4]
      end
      self.jacobGeoFunc = function (t, xn)
         return self.grid:calcJacobian(xn)
      end
   end

   self.bmagInvFunc = function (t, xn)
      return 1.0/self.bmagFunc(t,xn)
   end
   self.gradparFunc = function (t, xn)
      return 1.0/math.sqrt(self.g_zzFunc(t,xn))
   end
   self.geoXFunc = function (t, xn)
      return self.g_xzFunc(t,xn)/(self.bmagFunc(t,xn)*self.jacobGeoFunc(t,xn)*math.sqrt(self.g_zzFunc(t,xn)))
   end
   self.geoYFunc = function (t, xn)
      return self.g_yzFunc(t,xn)/(self.bmagFunc(t,xn)*self.jacobGeoFunc(t,xn)*math.sqrt(self.g_zzFunc(t,xn)))
   end
   self.geoZFunc = function (t, xn)
      return math.sqrt(self.g_zzFunc(t,xn))/(self.bmagFunc(t,xn)*self.jacobGeoFunc(t,xn))
   end

   -- use pseudo-local s-alpha geo
   if self.salpha then
      self.bmagInvFunc = function (t, xn)
         return 1.0/self.bmagFunc(t,xn)
      end
      self.gradparFunc = function (t, xn)
         return 1.0/self.salpha.q0/self.salpha.R0 
      end
      self.jacobGeoFunc = function (t, xn)
         return self.salpha.q0*self.salpha.R0
      end
      self.geoYFunc = function (t, xn)
         return -1.0/self.bmagFunc(t,xn)/self.salpha.r0
      end
      self.geoZFunc = function (t, xn)
         return 1.0/self.bmagFunc(t,xn)
      end
      self.gxx_Func = function (t, xn)
         return 1.0
      end
      if self.ndim == 1 then
         self.geoXFunc = function (t, xn)
            return self.salpha.shat*xn[1]/self.bmagFunc(t,xn)/self.salpha.r0
         end
         self.gxy_Func = function (t, xn)
            return self.salpha.shat*xn[1]
         end
         self.gyy_Func = function (t, xn)
            return 1.0 + self.salpha.shat^2*xn[1]^2
         end
      elseif self.ndim == 2 then 
         self.geoXFunc = function (t, xn)
            return 0.0
         end
         self.gxy_Func = function (t, xn)
            return 0.0
         end
         self.gyy_Func = function (t, xn)
            return 1.0
         end
      elseif self.ndim == 3 then
         self.geoXFunc = function (t, xn)
            return self.salpha.shat*xn[3]/self.bmagFunc(t,xn)/self.salpha.r0
         end
         self.gxy_Func = function (t, xn)
            return self.salpha.shat*xn[3]
         end
         self.gyy_Func = function (t, xn)
            return 1.0 + self.salpha.shat^2*xn[3]^2
         end
      end
   end

   -- inverse of jacobGeo, for removing jacobGeo factor from output quantities, e.g. density
   self.jacobGeoInvFunc = function (t, xn)
      return 1.0/self.jacobGeoFunc(t, xn)
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
   self.setGradpar = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gradparFunc
   }
   self.setJacobGeo = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobGeoFunc
   }
   self.setJacobGeoInv = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobGeoInvFunc
   }
   self.setGeoX = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoXFunc
   }
   self.setGeoY = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoYFunc
   }
   self.setGeoZ = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoZFunc
   }
   self.setGxx = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxx_Func
   }
   self.setGxy = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxy_Func
   }
   self.setGyy = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gyy_Func
   }
   if self.phiWallFunc then 
      self.setPhiWall = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         projectOnGhosts = true,
         evaluate = self.phiWallFunc
      }
   end

   -- determine which variables bmag depends on by checking if setting a variable to nan results in nan
   local gridCenter = {}
   for dir = 1, self.ndim do
      gridCenter[dir] = self.grid:mid(dir)
   end
   self.bmagVars = {}
   for dir = 1, self.ndim do
      gridCenter[dir] = 0/0 -- set this var to nan 
      -- test if result is nan.. nan is the only value that doesn't equal itself
      if self.bmagFunc(0, gridCenter) ~= self.bmagFunc(0, gridCenter) then 
        -- if result is nan, bmag must depend on this var
        table.insert(self.bmagVars, dir) 
      end
      gridCenter[dir] = self.grid:mid(dir) -- reset so we can check other vars
   end
   if self.bmagVars[1] == nil then self.bmagVars[1] = 0 end
end

function GkGeometry:createDiagnostics()
end

function GkGeometry:initField()
   self.setBmag:advance(0.0, {}, {self.geo.bmag})
   self.setBmagInv:advance(0.0, {}, {self.geo.bmagInv})
   self.setGradpar:advance(0.0, {}, {self.geo.gradpar})
   self.setJacobGeo:advance(0.0, {}, {self.geo.jacobGeo})
   self.setJacobGeoInv:advance(0.0, {}, {self.geo.jacobGeoInv})
   self.setGxx:advance(0.0, {}, {self.geo.gxx})
   self.setGxy:advance(0.0, {}, {self.geo.gxy})
   self.setGyy:advance(0.0, {}, {self.geo.gyy})
   self.setGeoX:advance(0.0, {}, {self.geo.geoX})
   self.setGeoY:advance(0.0, {}, {self.geo.geoY})
   self.setGeoZ:advance(0.0, {}, {self.geo.geoZ})
   if self.setPhiWall then self.setPhiWall:advance(0.0, {}, {self.geo.phiWall})
   else self.geo.phiWall:clear(0.0) end
   -- sync ghost cells. these calls do not enforce periodicity because
   -- these fields initialized with syncPeriodicDirs = false
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.gradpar:sync(false)
   self.geo.gxx:sync(false)
   self.geo.gxy:sync(false)
   self.geo.gyy:sync(false)
   self.geo.jacobGeo:sync(false)
   self.geo.jacobGeoInv:sync(false)
   self.geo.geoX:sync(false)
   self.geo.geoY:sync(false)
   self.geo.geoZ:sync(false)
   self.geo.phiWall:sync(false)

   -- apply BCs
   self:applyBc(0, self.geo)
end

function GkGeometry:write(tm)
   -- not evolving geometry, so only write geometry at beginning
   if self.ioFrame == 0 then
      self.fieldIo:write(self.geo.bmag, string.format("bmag_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.bmagInv, string.format("bmagInv_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gradpar, string.format("gradpar_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoX, string.format("geoX_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoY, string.format("geoY_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoZ, string.format("geoZ_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gxx, string.format("gxx_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gxy, string.format("gxy_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gyy, string.format("gyy_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.jacobGeo, string.format("jacobGeo_%d.bp", self.ioFrame), tm, self.ioFrame)
   end
   self.ioFrame = self.ioFrame+1
end

function GkGeometry:writeRestart(tm)
end

function GkGeometry:rkStepperFields()
   return { self.geo, self.geo, self.geo, self.geo }
end

function GkGeometry:advance(tCurr)
   if self.evolve then 
      self.setPhiWall:advance(tCurr, {}, self.geo.phiWall)
   end 
end

function GkGeometry:applyBc(tCurr, geoIn)
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
