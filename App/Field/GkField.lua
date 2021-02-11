-- GkField ---------------------------------------------------------------------
--
-- App support code: Gyrokinetic fields phi and apar, solved by
-- (perpendicular) Poisson and Ampere equations.
-- 
-- NOTE: GkGeometry is also in this file (farther down).
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Constants        = require "Lib.Constants"
local DataStruct       = require "DataStruct"
local LinearTrigger    = require "LinearTrigger"
local Mpi              = require "Comm.Mpi"
local Proto            = require "Lib.Proto"
local Updater          = require "Updater"
local xsys             = require "xsys"
local FieldBase        = require "App.Field.FieldBase"
local Species          = require "App.Species"
local Time             = require "Lib.Time"
local math             = require("sci.math").generic
local diff             = require("sci.diff")
local Logger           = require "Lib.Logger"
local lume             = require "Lib.lume"

local GkField = Proto(FieldBase.FieldBase)

local EM_BC_OPEN = 1

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkField:init(tbl)
   GkField.super.init(self, tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkField:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.
   
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- By default evolve field.

   self.isElectromagnetic = xsys.pickBool(tbl.isElectromagnetic, false) -- Electrostatic by default.

   -- Create triggers to write fields.
   local nFrame = tbl.nFrame or appTbl.nFrame
   self.ioTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   self.ioFrame = 0 -- Frame number for IO.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   
   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   -- Get boundary condition settings.
   -- These will be checked for consistency when the solver is initialized.
   if tbl.phiBcLeft then self.phiBcLeft          = tbl.phiBcLeft end
   if tbl.phiBcRight then self.phiBcRight        = tbl.phiBcRight end
   if tbl.phiBcBottom then self.phiBcBottom      = tbl.phiBcBottom end
   if tbl.phiBcTop then self.phiBcTop            = tbl.phiBcTop end
   if tbl.phiBcBack then self.phiBcBack          = tbl.phiBcBack end
   if tbl.phiBcFront then self.phiBcFront        = tbl.phiBcFront end
   if tbl.aparBcLeft then self.aparBcLeft        = tbl.aparBcLeft end
   if tbl.aparBcRight then self.aparBcRight      = tbl.aparBcRight end
   if tbl.aparBcBottom then self.aparBcBottom    = tbl.aparBcBottom end
   if tbl.aparBcTop then self.aparBcTop          = tbl.aparBcTop end
   if appTbl.periodicDirs then self.periodicDirs = appTbl.periodicDirs else self.periodicDirs = {} end

   -- For storing integrated energies.
   self.phiSq    = DataStruct.DynVector { numComponents = 1 }
   self.aparSq   = DataStruct.DynVector { numComponents = 1 }
   self.esEnergy = DataStruct.DynVector { numComponents = 1 }
   self.emEnergy = DataStruct.DynVector { numComponents = 1 }

   self.adiabatic = false
   self.discontinuousPhi  = xsys.pickBool(tbl.discontinuousPhi, false)
   self.discontinuousApar = xsys.pickBool(tbl.discontinuousApar, true)

   -- For ndim=1 only.
   self.kperp2 = tbl.kperp2

   -- Allow user to specify polarization weight. will be calculated automatically if not specified.
   self.polarizationWeight = tbl.polarizationWeight

   -- Determine whether to use linearized polarization term in poisson equation, which uses background density in polarization weight
   -- if not, uses full time-dependent density in polarization weight.
   self.linearizedPolarization = xsys.pickBool(tbl.linearizedPolarization, true)
   self.uniformPolarization    = xsys.pickBool(tbl.uniformPolarization, true)

   if self.isElectromagnetic then
      self.mu0 = tbl.mu0 or Constants.MU0
   end

   self.externalPhi = tbl.externalPhi
   if self.externalPhi and self.evolve then 
      print("GkField: warning... specifying externalPhi will make initial phi inconsistent with f") 
   end
   assert(not tbl.initPhiFunc, "GkField: initPhiFunc deprecated. Use externalPhi.")

   -- This allows us to apply a multiplicative time dependence to externalPhi.
   if tbl.externalPhiTimeDependence then
      self.externalPhiTimeDependence = tbl.externalPhiTimeDependence
   else
      self.externalPhiTimeDependence = false
   end

   -- Create trigger for how frequently to compute field energy.
   if appTbl.calcIntQuantEvery then
      self.calcIntFieldEnergyTrigger = LinearTrigger(0, appTbl.tEnd,  math.floor(1/appTbl.calcIntQuantEvery))
   else
      self.calcIntFieldEnergyTrigger = function(t) return true end
   end

   self.bcTime = 0.0 -- Timer for BCs.

   self._first     = true
   self._firstStep = true
end

-- Methods for EM field object.
function GkField:hasEB() return true, self.isElectromagnetic end
function GkField:setCfl() end
function GkField:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkField:setBasis(basis) self.basis = basis end
function GkField:setGrid(grid) self.grid = grid; self.ndim = self.grid:ndim() end

local function createField(grid, basis, ghostCells, vComp, periodicSync)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = ghostCells,
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      },
      syncPeriodicDirs = periodicSync,
   }
   fld:clear(0.0)
   return fld
end

function GkField:alloc(nRkDup)
   -- Allocate fields needed in RK update.
   -- nField is related to number of RK stages.
   self.potentials = {}
   self.nRkDup = nRkDup
   for i = 1, nRkDup do
      self.potentials[i] = {}
      self.potentials[i].phi    = createField(self.grid,self.basis,{1,1})
      self.potentials[i].phiAux = createField(self.grid,self.basis,{1,1})
      self.potentials[i].phi:clear(0.0)
      self.potentials[i].phiAux:clear(0.0)
      if self.isElectromagnetic then
         self.potentials[i].apar    = createField(self.grid,self.basis,{1,1})
         self.potentials[i].dApardt = createField(self.grid,self.basis,{1,1})
         self.potentials[i].apar:clear(0.0)
         self.potentials[i].dApardt:clear(0.0)
      end
   end

   self.dApardtProv = createField(self.grid,self.basis,{1,1})

   -- Create fields for total charge and current densities.
   self.chargeDens  = createField(self.grid,self.basis,{1,1})
   self.currentDens = createField(self.grid,self.basis,{1,1})
   -- Set up constant dummy field.
   self.unitWeight = createField(self.grid,self.basis,{1,1})
   local initUnit = Updater.ProjectOnBasis {
      onGrid   = self.grid,
      basis    = self.basis,
      evaluate = function (t,xn) return 1.0 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- Set up some other fields.
   self.weight          = createField(self.grid,self.basis,{1,1})
   self.laplacianWeight = createField(self.grid,self.basis,{1,1})
   self.modifierWeight  = createField(self.grid,self.basis,{1,1})
end

-- Solve for initial fields self-consistently 
-- from initial distribution function.
function GkField:initField(species)
   if self.externalPhi then
      local evalOnNodes = Updater.EvalOnNodes {
         onGrid          = self.grid,
         basis           = self.basis,
         evaluate        = self.externalPhi,
         projectOnGhosts = true
      }
      self.externalPhiFld = createField(self.grid,self.basis,{1,1})
      evalOnNodes:advance(0.0, {}, {self.externalPhiFld})
      for i = 1, self.nRkDup do
         if self.externalPhiTimeDependence then
            self.potentials[i].phi:combine(self.externalPhiTimeDependence(0.0),self.externalPhiFld)
         else
            self.potentials[i].phi:copy(self.externalPhiFld)
         end
      end
   else
      -- Solve for initial phi.
      self:advance(0.0, species, 1, 1)
   end

   local polyOrder = self.basis:polyOrder()

   if self.isElectromagnetic then
      -- Solve for initial Apar.
      local apar = self.potentials[1].apar
      self.currentDens:clear(0.0)
      for _, s in lume.orderedIter(species) do
         self.currentDens:accumulate(s:getCharge(), s:getMomDensity())
      end
      self.aparSlvr:advance(0.0, {self.currentDens}, {apar})

      -- Clear dApar/dt ... will be solved for before being used.
      self.potentials[1].dApardt:clear(0.0)
   end

   -- Apply BCs and update ghosts.
   self:applyBc(0, self.potentials[1])
end

function GkField:rkStepperFields()
   return self.potentials
end

-- For RK timestepping for non-elliptic fields (e.g. only apar).
function GkField:copyRk(outIdx, aIdx)
   if self.isElectromagnetic and self:rkStepperFields()[aIdx] then 
      self:rkStepperFields()[outIdx].apar:copy(self:rkStepperFields()[aIdx].apar)
   end
end
-- For RK timestepping for non-elliptic fields (e.g. only apar).
function GkField:combineRk(outIdx, a, aIdx, ...)
   if self.isElectromagnetic and self:rkStepperFields()[aIdx] then
      local args = {...} -- Package up rest of args as table.
      local nFlds = #args/2
      self:rkStepperFields()[outIdx].apar:combine(a, self:rkStepperFields()[aIdx].apar)
      for i = 1, nFlds do -- Accumulate rest of the fields.
         self:rkStepperFields()[outIdx].apar:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]].apar)
      end	 
   end
end

function GkField:suggestDt()
   return GKYL_MAX_DOUBLE
end
function GkField:clearCFL()
end

function GkField:createSolver(species, externalField)
   -- Get adiabatic species info.
   for _, s in lume.orderedIter(species) do
      if Species.AdiabaticSpecies.is(s) then
         self.adiabatic = true
         self.adiabSpec = s
      end
   end
   assert((self.adiabatic and self.isElectromagnetic) == false, "GkField: cannot use adiabatic response for electromagnetic case")

   -- Metric coefficients in the Poisson and Ampere equations for phi and Apar.
   local gxxPoisson, gxyPoisson, gyyPoisson, jacobGeo
   local gxxAmpere, gxyAmpere, gyyAmpere
   if externalField.geo then 
      if externalField.geo.name=="SimpleHelical" then
         -- Metric coefficients alone.
         gxxPoisson = externalField.geo.gxx
         gxyPoisson = externalField.geo.gxy
         gyyPoisson = externalField.geo.gyy
         gxxAmpere  = externalField.geo.gxx
         gxyAmpere  = externalField.geo.gxy
         gyyAmpere  = externalField.geo.gyy
         jacobGeo   = self.unitWeight
      elseif externalField.geo.name=="GenGeo" then
         -- Include Jacobian factor in metric coefficients (if using linearized polarization density).
         if self.linearizedPolarization then
            gxxPoisson = externalField.geo.gxxJ
            gxyPoisson = externalField.geo.gxyJ
            gyyPoisson = externalField.geo.gyyJ
         else  -- If not, Jacobian already included in polarization density or not needed.
            gxxPoisson = externalField.geo.gxx
            gxyPoisson = externalField.geo.gxy
            gyyPoisson = externalField.geo.gyy
         end
         gxxAmpere = externalField.geo.gxxJ
         gxyAmpere = externalField.geo.gxyJ
         gyyAmpere = externalField.geo.gyyJ
         jacobGeo  = externalField.geo.jacobGeo
      end
   else
      jacobGeo = self.unitWeight
   end
   self.phiSlvr = Updater.FemPoisson {
      onGrid   = self.grid,
      basis    = self.basis,
      bcLeft   = self.phiBcLeft,
      bcRight  = self.phiBcRight,
      bcBottom = self.phiBcBottom,
      bcTop    = self.phiBcTop,
      bcBack   = self.phiBcBack,
      bcFront  = self.phiBcFront,
      periodicDirs = self.periodicDirs,
      zContinuous  = not self.discontinuousPhi,
      -- Note these metric coefficients may contain the Jacobian (GenGeo).
      gxx = gxxPoisson,
      gxy = gxyPoisson,
      gyy = gyyPoisson,
   }
   if self.ndim == 3 and not self.discontinuousPhi then
      self.phiZSmoother = Updater.FemParPoisson {
        onGrid = self.grid,
        basis  = self.basis,
        smooth = true,
      }
   end
   -- when using a linearizedPolarization term in Poisson equation,
   -- The weights on the terms are constant scalars.
   if self.linearizedPolarization then
      -- If not provided, calculate species-dependent weight on polarization term == sum_s m_s n_s / B^2.
      if not self.polarizationWeight then 
         self.polarizationWeight = 0.0
         for _, s in lume.orderedIter(species) do
            if Species.GkSpecies.is(s) then
               self.polarizationWeight = self.polarizationWeight + s:getPolarizationWeight()
            end
         end
      end
      -- If not adiabatic, and polarization weight still not set, assume it is 1.
      if self.polarizationWeight == 0.0 and not self.adiabatic then self.polarizationWeight = 1.0 end

      -- Set up scalar multipliers on laplacian and modifier terms in Poisson equation.
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

      self.laplacianWeight:combine(laplacianConstant, self.unitWeight)   -- No jacobian here.
      self.modifierWeight:combine(modifierConstant, jacobGeo)

      if laplacianConstant ~= 0 then self.phiSlvr:setLaplacianWeight(self.laplacianWeight) end
      if modifierConstant ~= 0 then self.phiSlvr:setModifierWeight(self.modifierWeight) end
   end
   -- When not using linearizedPolarization, weights are set each step in advance method.

   if self.isElectromagnetic then 
     local ndim = self.ndim
     local laplacianConstant, modifierConstant
     -- NOTE: aparSlvr only used to solve for initial Apar
     -- at all other times Apar is timestepped using dApar/dt
     self.aparSlvr = Updater.FemPoisson {
        onGrid   = self.grid,
        basis    = self.basis,
        bcLeft   = self.aparBcLeft,
        bcRight  = self.aparBcRight,
        bcBottom = self.aparBcBottom,
        bcTop    = self.aparBcTop,
        periodicDirs = self.periodicDirs,
        zContinuous  = not self.discontinuousApar,
        -- Note these metric coefficients may contain the Jacobian (GenGeo).
        gxx = gxxAmpere,
        gxy = gxyAmpere,
        gyy = gyyAmpere,
     }
     if ndim==1 then
        laplacianConstant = 0.0
        modifierConstant  = self.kperp2/self.mu0
     else
        laplacianConstant = -1.0/self.mu0
        modifierConstant  = 0.0
     end
     self.laplacianWeight:combine(laplacianConstant, self.unitWeight)
     self.modifierWeight:combine(modifierConstant, self.unitWeight)
     if laplacianConstant ~= 0 then self.aparSlvr:setLaplacianWeight(self.laplacianWeight) end
     if modifierConstant ~= 0 then self.aparSlvr:setModifierWeight(self.modifierWeight) end

     self.dApardtSlvr = Updater.FemPoisson {
        onGrid   = self.grid,
        basis    = self.basis,
        bcLeft   = self.aparBcLeft,
        bcRight  = self.aparBcRight,
        bcBottom = self.aparBcBottom,
        bcTop    = self.aparBcTop,
        periodicDirs = self.periodicDirs,
        zContinuous  = not self.discontinuousApar,
        -- Note these metric coefficients may contain the Jacobian (GenGeo).
        gxx = gxxAmpere,
        gxy = gxyAmpere,
        gyy = gyyAmpere,
     }
     if ndim==1 then
        laplacianConstant = 0.0
        modifierConstant  = 1.0
     else
        laplacianConstant = -1.0/self.mu0
        modifierConstant  = 1.0
     end
     self.laplacianWeight:combine(laplacianConstant, self.unitWeight)
     self.modifierWeight:combine(modifierConstant, self.unitWeight)
     if laplacianConstant ~= 0 then self.dApardtSlvr:setLaplacianWeight(self.laplacianWeight) end
     if modifierConstant ~= 0 then self.dApardtSlvr:setModifierWeight(self.modifierWeight) end

     -- Separate solver for additional step for p=1.
     if self.basis:polyOrder() == 1 then
        self.dApardtSlvr2 = Updater.FemPoisson {
           onGrid   = self.grid,
           basis    = self.basis,
           bcLeft   = self.aparBcLeft,
           bcRight  = self.aparBcRight,
           bcBottom = self.aparBcBottom,
           bcTop    = self.aparBcTop,
           periodicDirs = self.periodicDirs,
           zContinuous  = not self.discontinuousApar,
           -- Note these metric coefficients may contain the Jacobian (GenGeo).
           gxx = gxxAmpere,
           gxy = gxyAmpere,
           gyy = gyyAmpere,
        }
        if ndim==1 then
           laplacianConstant = 0.0
           modifierConstant  = 1.0
        else
           laplacianConstant = -1.0/self.mu0
           modifierConstant  = 1.0
        end
        self.laplacianWeight:combine(laplacianConstant, self.unitWeight)
        self.modifierWeight:combine(modifierConstant, self.unitWeight)
        if laplacianConstant ~= 0 then self.dApardtSlvr2:setLaplacianWeight(self.laplacianWeight) end
        if modifierConstant ~= 0 then self.dApardtSlvr2:setModifierWeight(self.modifierWeight) end
     end
   end

   -- Need to set this flag so that field calculated self-consistently at end of full RK timestep.
   self.isElliptic = true

   if self.isElectromagnetic and self.basis:polyOrder() == 1 then 
      self.nstep = 3 
   elseif self.isElectromagnetic then
      self.nstep = 2
   end

   -- Function to construct a BC updater.
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
	 onGrid = self.grid,
	 boundaryConditions = bcList,
	 dir = dir,
	 edge = edge,
      }
   end

   -- Various functions to apply BCs of different types.
   local function bcOpen(dir, tm, idxIn, fIn, fOut)
      -- Requires skinLoop = "pointwise".
      self.basis:flipSign(dir, fIn, fOut)
   end

   -- Functions to make life easier while reading in BCs to apply.
   self.boundaryConditions = { }   -- List of Bcs to apply.
   local function appendBoundaryConditions(dir, edge, bcType)
      if bcType == EM_BC_OPEN then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, { bcOpen }))
      else
	 assert(false, "GkField: Unsupported BC type!")
      end
   end

   local function handleBc(dir, bc)
      if bc[1] then
	 appendBoundaryConditions(dir, "lower", bc[1])
      end
      if bc[2] then
	 appendBoundaryConditions(dir, "upper", bc[2])
      end
   end

   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

   -- For non-periodic dirs, use BC_OPEN to make sure values on edge of ghost cells match
   -- values on edge of skin cells, so that field is continuous across skin-ghost boundary.
   for dir = 1, self.ndim do
      if not contains(self.periodicDirs, dir) then 
         handleBc(dir, {EM_BC_OPEN, EM_BC_OPEN}) 
      end
   end
end

function GkField:createDiagnostics()
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.potentials[1].phi:elemType(),
      method     = self.ioMethod,
      writeGhost = self.writeGhost,
      metaData   = {
	 polyOrder = self.basis:polyOrder(),
	 basisType = self.basis:id()
      },
   }

   -- Updaters for computing integrated quantities.
   self.int2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid   = self.grid,
      basis    = self.basis,
      quantity = "V2"
   }
   if self.ndim == 1 then
      self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
         onGrid   = self.grid,
         basis    = self.basis,
         quantity = "V2"
      }
   elseif self.basis:polyOrder()==1 then -- GradPerpV2 only implemented for p=1 currently.
      self.energyCalc = Updater.CartFieldIntegratedQuantCalc {
         onGrid   = self.grid,
         basis    = self.basis,
         quantity = "GradPerpV2"
      }
   end
end

function GkField:write(tm, force)
   if self.evolve then
      if self.calcIntFieldEnergyTrigger(tm) then
         -- Compute integrated quantities over domain.
         self.int2Calc:advance(tm, { self.potentials[1].phi }, { self.phiSq })
         if self.isElectromagnetic then 
            self.int2Calc:advance(tm, { self.potentials[1].apar }, { self.aparSq })
         end
         if self.energyCalc then 
            if self.linearizedPolarization then
               local esEnergyFac = .5*self.polarizationWeight
               if self.ndim == 1 then 
                  esEnergyFac = esEnergyFac*self.kperp2 
                  if self.adiabatic then 
                     esEnergyFac = esEnergyFac + .5*self.adiabSpec:getQneutFac() 
                  end
               end
               self.energyCalc:advance(tm, { self.potentials[1].phiAux, esEnergyFac }, { self.esEnergy })

               -- NRM: the section below adds a (physical) term to esEnergy proportional to |phi|**2 when using adiabatic electrons.
               -- however, this makes it difficult to compute the growth rate of ky!=0 modes from esEnergy when the
               -- ky=0 mode is not suppressed. commenting out the below is a hack so that esEnergy can be used to compute
               -- growth rate of ky!=0 modes.

               --if self.adiabatic and self.ndim > 1 then
               --   local tm, energyVal = self.esEnergy:lastData()
               --   local _, phiSqVal = self.phiSq:lastData()
               --   energyVal[1] = energyVal[1] + .5*self.adiabSpec:getQneutFac()*phiSqVal[1]
               --end
            else
               -- Something.
            end
            if self.isElectromagnetic then 
              local emEnergyFac = .5/self.mu0
              if self.ndim == 1 then emEnergyFac = emEnergyFac*self.kperp2 end
              self.energyCalc:advance(tm, { self.potentials[1].apar, emEnergyFac}, { self.emEnergy })
            end
         end
      end

      if self.ioFrame == 0 then 
         self.fieldIo:write(self.phiSlvr:getLaplacianWeight(), "laplacianWeight_0.bp", tm, self.ioFrame, false)
         self.fieldIo:write(self.phiSlvr:getModifierWeight(), "modifierWeight_0.bp", tm, self.ioFrame, false)
      end
      
      if self.ioTrigger(tm) or force then
	 self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
         if self.isElectromagnetic then 
	   self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
	   self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
	 self.phiSq:write(string.format("phiSq.bp"), tm, self.ioFrame)
	 self.esEnergy:write(string.format("esEnergy.bp"), tm, self.ioFrame)
	 if self.isElectromagnetic then
	    self.aparSq:write(string.format("aparSq.bp"), tm, self.ioFrame)
	    self.emEnergy:write(string.format("emEnergy.bp"), tm, self.ioFrame)
	 end
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
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
   -- (the final "false" prevents writing of ghost cells).
   self.fieldIo:write(self.potentials[1].phi, "phi_restart.bp", tm, self.ioFrame, false)
   self.fieldIo:write(self.phiSlvr:getLaplacianWeight(), "laplacianWeight_restart.bp", tm, self.ioFrame, false)
   self.fieldIo:write(self.phiSlvr:getModifierWeight(), "modifierWeight_restart.bp", tm, self.ioFrame, false)
   if self.isElectromagnetic then
     self.fieldIo:write(self.potentials[1].apar, "apar_restart.bp", tm, self.ioFrame, false)
   end

   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self.phiSq:write("phiSq_restart.bp", tm, self.dynVecRestartFrame, false, false)
   self.dynVecRestartFrame = self.dynVecRestartFrame + 1

end

function GkField:readRestart()
   -- This read of restart file of phi is only to get frame
   -- numbering correct. The forward Euler recomputes the potential
   -- before updating the hyperbolic part.
   local tm, fr = self.fieldIo:read(self.potentials[1].phi, "phi_restart.bp")
   self.phiSq:read("phiSq_restart.bp", tm)

   self.fieldIo:read(self.laplacianWeight, "laplacianWeight_restart.bp")
   self.fieldIo:read(self.modifierWeight, "modifierWeight_restart.bp")
   self.phiSlvr:setLaplacianWeight(self.laplacianWeight)
   if self.adiabatic or self.ndim == 1 then self.phiSlvr:setModifierWeight(self.modifierWeight) end

   if self.isElectromagnetic then
      self.fieldIo:read(self.potentials[1].apar, "apar_restart.bp")
   end

   self:applyBc(0, self.potentials[1])

   self.ioFrame = fr 
   -- Iterate triggers.
   self.ioTrigger(tm)
end

-- Not needed for GK.
function GkField:accumulateCurrent(dt, current, em)
end

-- Solve for electrostatic potential phi.
function GkField:advance(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs  = self:rkStepperFields()[outIdx]
   
   if self.evolve or (self._first and not self.externalPhi) then
      if self.externalPhiTimeDependence then
         potCurr.phi:combine(self.externalPhiTimeDependence(tCurr), self.externalPhiFld)
      else
         self.chargeDens:clear(0.0)
         for _, s in lume.orderedIter(species) do
            self.chargeDens:accumulate(s:getCharge(), s:getNumDensity())
         end
         -- If not using linearized polarization term, set up laplacian weight.
         if not self.linearizedPolarization or (self._first and not self.uniformPolarization) then
            self.weight:clear(0.0)
            for _, s in lume.orderedIter(species) do
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
         -- Phi solve (elliptic, so update potCurr).
         -- Energy conservation requires phi to be continuous in all directions. 
         -- The first FEM solve ensures that phi is continuous in x and y.
         -- The conserved energy is defined in terms of this intermediate result,
         -- which we denote phiAux, before the final smoothing operation in z.
         self.phiSlvr:advance(tCurr, {self.chargeDens}, {potCurr.phiAux})
         -- Smooth phi in z to ensure continuity in all directions.
         if self.ndim == 3 and not self.discontinuousPhi then
            self.phiZSmoother:advance(tCurr, {potCurr.phiAux}, {potCurr.phi})
         else
            potCurr.phi = potCurr.phiAux
         end

         -- Apply BCs.
         local tmStart = Time.clock()
         -- make sure phi is continuous across skin-ghost boundary
         for _, bc in ipairs(self.boundaryConditions) do
            bc:advance(tCurr, {}, {potCurr.phi})
         end
         potCurr.phi:sync(true)
         self.bcTime = self.bcTime + (Time.clock()-tmStart)
 
         self._first = false
      end
   else
      -- Just copy stuff over.
      if self.isElectromagnetic then 
         potRhs.apar:copy(potCurr.apar) 
      end
   end
end

-- Solve for dApardt in p>=2, or solve for a provisional dApardtProv in p=1.
function GkField:advanceStep2(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs  = self:rkStepperFields()[outIdx]

   local polyOrder = self.basis:polyOrder()

   if self.evolve then

      self.currentDens:clear(0.0)
      if self.ndim==1 then 
         self.modifierWeight:combine(self.kperp2/self.mu0, self.unitWeight) 
      else 
         self.modifierWeight:clear(0.0)
      end
      for _, s in lume.orderedIter(species) do
         if s:isEvolving() then 
            self.modifierWeight:accumulate(s:getCharge()*s:getCharge()/s:getMass(), s:getNumDensity())
            -- Taking momDensity at outIdx gives momentum moment of df/dt.
            self.currentDens:accumulate(s:getCharge(), s:getMomDensity(outIdx))
         end
      end
      self.dApardtSlvr:setModifierWeight(self.modifierWeight)
      -- dApar/dt solve.
      local dApardt = potCurr.dApardt
      self.dApardtSlvr:advance(tCurr, {self.currentDens}, {dApardt}) 

      -- Apply BCs.
      local tmStart = Time.clock()
      -- make sure dApardt is continuous across skin-ghost boundary
      for _, bc in ipairs(self.boundaryConditions) do
         bc:advance(tCurr, {}, {potCurr.dApardt})
      end
      dApardt:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      if polyOrder > 1 then 
         -- Apar RHS is just dApar/dt.
         potRhs.apar:copy(dApardt)
      else 
         -- Save dApardt as dApardtProv, so that it can be used in upwinding 
         -- in vpar surface terms in p=1 Ohm's law and GK update.
         self.dApardtProv:copy(dApardt)
      end
   end
end

-- Note: step 3 is for p=1 only: solve for dApardt.
function GkField:advanceStep3(tCurr, species, inIdx, outIdx)
   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs  = self:rkStepperFields()[outIdx]

   local polyOrder = self.basis:polyOrder()

   if self.evolve then
      self.currentDens:clear(0.0)
      if self.ndim==1 then 
         self.modifierWeight:combine(self.kperp2/self.mu0, self.unitWeight) 
      else 
         self.modifierWeight:clear(0.0)
      end
      for _, s in lume.orderedIter(species) do
         if s:isEvolving() then 
            self.modifierWeight:accumulate(s:getCharge()*s:getCharge()/s:getMass(), s:getEmModifier())
            -- Taking momDensity at outIdx gives momentum moment of df/dt.
            self.currentDens:accumulate(s:getCharge(), s:getMomProjDensity(outIdx))
         end
      end
      self.dApardtSlvr2:setModifierWeight(self.modifierWeight)
      -- dApar/dt solve.
      local dApardt = potCurr.dApardt
      self.dApardtSlvr2:advance(tCurr, {self.currentDens}, {dApardt}) 

      -- Apply BCs.
      local tmStart = Time.clock()
      -- make sure dApardt is continuous across skin-ghost boundary
      for _, bc in ipairs(self.boundaryConditions) do
         bc:advance(tCurr, {}, {potCurr.dApardt})
      end
      dApardt:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      -- Apar RHS is just dApar/dt.
      potRhs.apar:copy(dApardt)
   end
end

function GkField:applyBcIdx(tCurr, idx)
   -- don't do anything here. BCs handled in advance steps.
end

-- NOTE: global boundary conditions handled by solver. this just updates interproc ghosts.
-- Also NOTE: this method does not usually get called (because it is not called in applyBcIdx).
function GkField:applyBc(tCurr, potIn)
   local tmStart = Time.clock()
   for _, bc in ipairs(self.boundaryConditions) do
      bc:advance(tCurr, {}, {potIn.phi})
   end
   potIn.phi:sync(true)
   if self.isElectromagnetic then 
     -- make sure apar is continuous across skin-ghost boundary
     for _, bc in ipairs(self.boundaryConditions) do
        bc:advance(tCurr, {}, {potIn.apar})
     end
     potIn.apar:sync(true) 
     -- make sure dApardt is continuous across skin-ghost boundary
     for _, bc in ipairs(self.boundaryConditions) do
        bc:advance(tCurr, {}, {potIn.dApardt})
     end
     potIn.dApardt:sync(true) 
   end
   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end
   
function GkField:totalSolverTime()
   if self.phiSlvr then
     local time = self.phiSlvr.totalTime
     if self.isElectromagnetic and self.aparSlvr then time = time + self.aparSlvr.totalTime end
     if self.isElectromagnetic and self.dApardtSlvr then time = time + self.dApardtSlvr.totalTime end
     if self.isElectromagnetic and self.dApardtSlvr2 then time = time + self.dApardtSlvr2.totalTime end
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

function GkField:printDevDiagnostics()
   self.phiSlvr:printDevDiagnostics()
end

-- GkGeometry ---------------------------------------------------------------------
--
-- A field object with fields specifying the magnetic geometry for GK.
--
--------------------------------------------------------------------------------

local GkGeometry = Proto(FieldBase.ExternalFieldBase)

-- Methods for no field object.
function GkGeometry:init(tbl)
   GkGeometry.super.init(self, tbl)
   self.tbl = tbl
end

function GkGeometry:fullInit(appTbl)
   local tbl = self.tbl -- previously store table.

   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, false) -- by default these fields are not time-dependent.

   -- Create triggers to write fields.
   local nFrame = tbl.nFrame or appTbl.nFrame
   self.ioTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   self.ioFrame = 0 -- Frame number for IO.
   
   -- Get function to initialize background magnetic field.
   self.bmagFunc = tbl.bmag
   --assert(self.bmagFunc and type(self.bmagFunc)=="function", "GkGeometry: must specify background magnetic field function with 'bmag'")

   -- Specify which geometry to use. This gets redefined in :alloc method if mapc2p is detected.
   -- NOTE: for must purposes one must use geo.name, defined in the :alloc method.
   self.geoType = tbl.geometryType and tbl.geometryType or "SimpleHelical"

   -- Specify which geometry to use.
   self.geoType = tbl.geometryType and tbl.geometryType or "SimpleHelical"

   -- Wall potential for sheath BCs.
   self.phiWallFunc = tbl.phiWall
   if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "GkGeometry: phiWall must be a function (t, xn)") end

   -- File containing geometry quantities that go into equations.
   self.fromFile = tbl.fromFile

   -- write ghost cells on boundaries of global domain (for BCs)
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)
end

function GkGeometry:hasEB() end
function GkGeometry:setCfl() end
function GkGeometry:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function GkGeometry:setBasis(basis) self.basis = basis end
function GkGeometry:setGrid(grid) self.grid = grid; self.ndim = self.grid:ndim() end

function GkGeometry:alloc()
   -- Allocate fields.
   self.geo = {}

   -- If mapc2p is detected, use general geometry infrastructure.
   if self.grid._mapc2p then
      self.geo.name = "GenGeo"
   else
      self.geo.name = self.geoType
   end

   local ghostNum     = {1,1}
   local syncPeriodic = false

   -- Background magnetic field.
   self.geo.bmag = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- bmagInv ~ 1/B.
   self.geo.bmagInv = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- gradpar = J B / sqrt(g_zz).
   self.geo.gradpar = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- Functions for Laplacian.
   -- g^xx = |grad x|**2.
   self.geo.gxx = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- g^xy = grad x . grad y.
   self.geo.gxy = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- g^yy = |grad y|**2.
   self.geo.gyy = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   if self.geo.name == "SimpleHelical" then

      -- Functions for magnetic drifts.
      -- bdriftX = 1/B*curl(bhat).grad x.
      self.geo.bdriftX = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
      -- bdriftY = 1/B*curl(bhat).grad y.
      self.geo.bdriftY = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      if self.ndim == 1 then
         self.geo.allGeo = createField(self.grid,self.basis,ghostNum,6,syncPeriodic)
      else
         self.geo.allGeo = createField(self.grid,self.basis,ghostNum,8,syncPeriodic)
      end

   elseif self.geo.name == "GenGeo" then

      -- Jacobian of coordinate transformation.
      self.geo.jacobGeo = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      -- Inverse of jacobian of coordinate transformation.
      self.geo.jacobGeoInv = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      -- Total jacobian, including phase space jacobian.
      -- jacobTot = J B 
      self.geo.jacobTot = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      -- Inverse of total jacobian, including phase space jacobian.
      -- jacobTotInv = 1 / ( J B )
      self.geo.jacobTotInv = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      -- Functions for magnetic drifts .
      -- geoX = g_xz / ( sqrt(g_zz) )
      self.geo.geoX = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
      -- geoY = g_yz / ( sqrt(g_zz) )
      self.geo.geoY = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
      -- geoZ = sqrt(g_zz)
      self.geo.geoZ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
 
      -- Functions for laplacian, including Jacobian factor.
      self.geo.gxxJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
      self.geo.gxyJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
      self.geo.gyyJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   
      if self.fromFile == nil then
         self.geo.allGeo = createField(self.grid,self.basis,ghostNum,16,syncPeriodic)
      end

   end

   -- Wall potential for sheath BCs.
   self.geo.phiWall = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

      
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.geo.bmag:elemType(),
      method     = self.ioMethod,
      writeGhost = self.writeGhost,
      metaData   = {
	 polyOrder = self.basis:polyOrder(),
	 basisType = self.basis:id()
      },
   }   
end

function GkGeometry:createSolver()

   -- Get configuration-space Jacobian from grid (used by App/Projection/GkProjection, via GkSpecies).
   self.jacobGeoFunc = function (t, xn)
      return self.grid:calcJacobian(xn)
   end

   if self.geo.name == "SimpleHelical" then

      -- Calculate magnetic drift functions.
      if self.ndim > 1 then
         local function bgrad(xn)
            local function bmagUnpack(...)
               local xn1 = {...}
               return self.bmagFunc(0, xn1)
            end
            local deriv = diff.derivativef(bmagUnpack, #xn)
            local xntable = {}
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

      -- Calculate all geometry quantities at once.
      if self.ndim == 1 then
         self.calcAllGeo = function(t, xn)
            local gxx, gxy, gyy = 1.0, 0.0, 1.0

            local bmag    = self.bmagFunc(t, xn)
            local gradpar = 1.0

            return bmag, 1/bmag, gradpar, gxx, gxy, gyy
          end
      else
         self.calcAllGeo = function(t, xn)
            local gxx, gxy, gyy = 1.0, 0.0, 1.0

            local bmag    = self.bmagFunc(t, xn)
            local gradpar = 1.0

            local bdriftX = self.bdriftXFunc(t, xn)
            local bdriftY = self.bdriftYFunc(t, xn)

            return bmag, 1/bmag, gradpar, gxx, gxy, gyy, bdriftX, bdriftY
          end
      end

      if self.fromFile == nil then
         self.setAllGeo = Updater.ProjectOnBasis {
            onGrid   = self.grid,
            basis    = self.basis,
            evaluate = self.calcAllGeo,
            projectOnGhosts = true,
         }
      end

      -- Determine which variables bmag depends on by checking if setting a variable to nan results in nan.
      local ones = {}
      for dir = 1, self.ndim do ones[dir] = 1 end
      self.bmagVars = {}
      for dir = 1, self.ndim do
         ones[dir] = 0/0 -- Set this var to nan.
         -- Test if result is nan.. nan is the only value that doesn't equal itself.
         if self.bmagFunc(0, ones) ~= self.bmagFunc(0, ones) then
           -- If result is nan, bmag must depend on this var.
           table.insert(self.bmagVars, dir)
         end
         ones[dir] = 1 -- Reset so we can check other vars.
      end
      if self.bmagVars[1] == nil then self.bmagVars[1] = 0 end

   elseif self.geo.name == "GenGeo" then

      -- Calculate all geometry quantities at once to avoid repeated metric calculations.
      if self.ndim == 1 then
         self.calcAllGeo = function(t, xn)
            local g = {}
            self.grid:calcMetric(xn, g)
            local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
            if self.grid._inDim==1 then
               g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = 1.0, 0.0, 0.0, 1.0, 0.0, g[1]
            elseif self.grid._inDim==2 then
               g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], 0.0, g[3], 0.0, 1.0
            elseif self.grid._inDim==3 then
               g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], g[3], g[4], g[5], g[6]
            end
            local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
            local geoX     = g_xz/math.sqrt(g_zz)
            local geoY     = g_yz/math.sqrt(g_zz)
            local geoZ     = math.sqrt(g_zz)

            local det = jacobian^2
            local gxx = (g_yy*g_zz-g_yz^2)/det
            local gxy = (g_xz*g_yz-g_xy*g_zz)/det
            local gxz = (g_xy*g_yz-g_xz*g_yy)/det
            local gyy = (g_xx*g_zz-g_xz^2)/det
            local gyz = (g_xy*g_xz-g_xx*g_yz)/det
            local gzz = (g_xx*g_yy-g_xy^2)/det

            local bmag    = self.bmagFunc(t, xn)
            local gradpar = jacobian*bmag/math.sqrt(g_zz)

            return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, 1/bmag, gradpar, 
                   geoX, geoY, geoZ, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
         end
      elseif self.ndim == 2 then
         self.calcAllGeo = function(t, xn)
            local g = {}
            self.grid:calcMetric(xn, g)
            local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
            if self.grid._inDim==2 then
               g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], 0.0, g[3], 0.0, 1.0
            elseif self.grid._inDim==3 then
               g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], g[3], g[4], g[5], g[6]
            end
            local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
            local geoX     = g_xz/math.sqrt(g_zz)
            local geoY     = g_yz/math.sqrt(g_zz)
            local geoZ     = math.sqrt(g_zz)

            local det = jacobian^2
            local gxx = (g_yy*g_zz-g_yz^2)/det
            local gxy = (g_xz*g_yz-g_xy*g_zz)/det
            local gxz = (g_xy*g_yz-g_xz*g_yy)/det
            local gyy = (g_xx*g_zz-g_xz^2)/det
            local gyz = (g_xy*g_xz-g_xx*g_yz)/det
            local gzz = (g_xx*g_yy-g_xy^2)/det

            local bmag    = self.bmagFunc(t, xn)
            local gradpar = jacobian*bmag/math.sqrt(g_zz)

            return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, 1/bmag, gradpar, 
                   geoX, geoY, geoZ, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
          end
      else
         self.calcAllGeo = function(t, xn)
            local g = {}
            self.grid:calcMetric(xn, g)
            local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], g[3], g[4], g[5], g[6]
            local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
            local geoX     = g_xz/math.sqrt(g_zz)
            local geoY     = g_yz/math.sqrt(g_zz)
            local geoZ     = math.sqrt(g_zz)

            local det = jacobian^2
            local gxx = (g_yy*g_zz-g_yz^2)/det
            local gxy = (g_xz*g_yz-g_xy*g_zz)/det
            local gxz = (g_xy*g_yz-g_xz*g_yy)/det
            local gyy = (g_xx*g_zz-g_xz^2)/det
            local gyz = (g_xy*g_xz-g_xx*g_yz)/det
            local gzz = (g_xx*g_yy-g_xy^2)/det

            local bmag    = self.bmagFunc(t, xn)
            local gradpar = jacobian*bmag/math.sqrt(g_zz)

            return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, 1/bmag, gradpar, 
                   geoX, geoY, geoZ, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
          end
      end

      if self.ndim == 3 then
         self.bmagVars = {1,3} 
      else
         self.bmagVars = {1}
      end

      if self.fromFile == nil then
         self.setAllGeo = Updater.EvalOnNodes {
            onGrid   = self.grid,
            basis    = self.basis,
            evaluate = self.calcAllGeo,
            projectOnGhosts = true,
         }
      end

   end   -- End of if geo.name == "..." statement.

   if self.phiWallFunc then 
      self.setPhiWall = Updater.EvalOnNodes {
         onGrid          = self.grid,
         basis           = self.basis,
         projectOnGhosts = true,
         evaluate        = self.phiWallFunc
      }
   end


   self.unitWeight = createField(self.grid,self.basis,{1,1})
   local initUnit = Updater.ProjectOnBasis {
      onGrid   = self.grid,
      basis    = self.basis,
      evaluate = function (t,xn) return 1.0 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- Updater to separate vector components packed into a single CartField.
   self.separateComponents = Updater.SeparateVectorComponents {
      onGrid = self.grid,
      basis  = self.basis,
   }
end

function GkGeometry:createDiagnostics()
end

function GkGeometry:initField()
   if self.geo.name == "SimpleHelical" then
      if self.fromFile then
         -- Read the geometry quantities from a file.
         if self.ndim == 1 then
            local tm, fr = self.fieldIo:read({bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
               gradpar=self.geo.gradpar, gxx=self.geo.gxx, gxy=self.geo.gxy, gyy=self.geo.gyy}, self.fromFile)
         else
            local tm, fr = self.fieldIo:read({bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
               gradpar=self.geo.gradpar, gxx=self.geo.gxx, gxy=self.geo.gxy, gyy=self.geo.gyy, 
               bdriftX=self.geo.bdriftX, bdriftY=self.geo.bdriftY}, self.fromFile)
         end
      else
         self.setAllGeo:advance(0.0, {}, {self.geo.allGeo})
         if self.ndim == 1 then
            self.separateComponents:advance(0, {self.geo.allGeo},
               {self.geo.bmag, self.geo.bmagInv, self.geo.gradpar,
                self.geo.gxx, self.geo.gxy, self.geo.gyy})
         else
            self.separateComponents:advance(0, {self.geo.allGeo},
               {self.geo.bmag, self.geo.bmagInv, self.geo.gradpar,
                self.geo.gxx, self.geo.gxy, self.geo.gyy, self.geo.bdriftX, self.geo.bdriftY})
         end
      end
   elseif self.geo.name == "GenGeo" then
      if self.fromFile then
         -- Read the geometry quantities from a file.
         local tm, fr = self.fieldIo:read({jacobGeo=self.geo.jacobGeo, jacobGeoInv=self.geo.jacobGeoInv, jacobTot=self.geo.jacobTot,
            jacobTotInv=self.geo.jacobTotInv, bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
            gradpar=self.geo.gradpar, geoX=self.geo.geoX, geoY=self.geo.geoY, geoZ=self.geo.geoZ, gxx=self.geo.gxx,
            gxy=self.geo.gxy, gyy=self.geo.gyy, gxxJ=self.geo.gxxJ, gxyJ=self.geo.gxyJ, gyyJ=self.geo.gyyJ},
            self.fromFile)
      else
         self.setAllGeo:advance(0.0, {}, {self.geo.allGeo})
         self.separateComponents:advance(0, {self.geo.allGeo},
            {self.geo.jacobGeo, self.geo.jacobGeoInv, self.geo.jacobTot, self.geo.jacobTotInv,
             self.geo.bmag, self.geo.bmagInv, self.geo.gradpar, self.geo.geoX, self.geo.geoY, self.geo.geoZ,
             self.geo.gxx, self.geo.gxy, self.geo.gyy, self.geo.gxxJ, self.geo.gxyJ, self.geo.gyyJ})
      end
   end
   local log = Logger { logToFile = true }
   log("Finished initializing the geometry\n")

   if self.setPhiWall then self.setPhiWall:advance(0.0, {}, {self.geo.phiWall})
   else self.geo.phiWall:clear(0.0) end

   -- Sync ghost cells. These calls do not enforce periodicity because
   -- these fields were created with syncPeriodicDirs = false.
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.gradpar:sync(false)
   self.geo.gxx:sync(false)
   self.geo.gxy:sync(false)
   self.geo.gyy:sync(false)
   if self.geo.name == "SimpleHelical" then
      if self.ndim > 1 then
         self.geo.bdriftX:sync(false)
         self.geo.bdriftY:sync(false)
      end
   elseif self.geo.name == "GenGeo" then
      self.geo.gxxJ:sync(false)
      self.geo.gxyJ:sync(false)
      self.geo.gyyJ:sync(false)
      self.geo.jacobGeo:sync(false)
      self.geo.jacobGeoInv:sync(false)
      self.geo.jacobTotInv:sync(false)
      self.geo.geoX:sync(false)
      self.geo.geoY:sync(false)
      self.geo.geoZ:sync(false)
   end
   self.geo.phiWall:sync(false)

   -- Apply BCs.
   self:applyBc(0, self.geo)
end

function GkGeometry:write(tm)
   -- Not evolving geometry, so only write geometry at beginning.
   if self.ioFrame == 0 then
      -- Write the geometry quantities to a file.
      if self.geo.name == "SimpleHelical" then
         if self.ndim == 1 then
            self.fieldIo:write({bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
               gradpar=self.geo.gradpar, gxx=self.geo.gxx, gxy=self.geo.gxy, gyy=self.geo.gyy},
               string.format("allGeo_%d.bp", self.ioFrame), tm, self.ioFrame)
         else
            self.fieldIo:write({bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
               gradpar=self.geo.gradpar, gxx=self.geo.gxx, gxy=self.geo.gxy, gyy=self.geo.gyy,
               bdriftX=self.geo.bdriftX, bdriftY=self.geo.bdriftY},
               string.format("allGeo_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
      elseif self.geo.name == "GenGeo" then
         self.fieldIo:write({jacobGeo=self.geo.jacobGeo, jacobGeoInv=self.geo.jacobGeoInv, jacobTot=self.geo.jacobTot,
            jacobTotInv=self.geo.jacobTotInv, bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
            gradpar=self.geo.gradpar, geoX=self.geo.geoX, geoY=self.geo.geoY, geoZ=self.geo.geoZ, gxx=self.geo.gxx,
            gxy=self.geo.gxy, gyy=self.geo.gyy, gxxJ=self.geo.gxxJ, gxyJ=self.geo.gxyJ, gyyJ=self.geo.gyyJ},
            string.format("allGeo_%d.bp", self.ioFrame), tm, self.ioFrame)

         -- Write a grid file.
         local metaData = {
            polyOrder = self.basis:polyOrder(),
            basisType = self.basis:id(),
            grid      = GKYL_OUT_PREFIX .. "_grid.bp"
         }
         self.grid:write("grid.bp", 0.0, metaData)
      end

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

function GkGeometry:applyBcIdx(tCurr, idx)
   self:applyBc(tCurr, self:rkStepperFields()[1])
end

function GkGeometry:applyBc(tCurr, geoIn)
   if self.evolve then geoIn.phiWall:sync(false) end
end

function GkGeometry:totalSolverTime()
   if self.evolve then return self.setPhiWall.totalTime
   else return 0 end
end

function GkGeometry:totalBcTime() return 0.0 end
function GkGeometry:energyCalcTime() return 0.0 end


return {GkField = GkField, GkGeometry = GkGeometry}
