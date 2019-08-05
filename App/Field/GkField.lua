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
local Mpi = require "Comm.Mpi"
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

local EM_BC_OPEN = 1

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

   local polyOrder = self.basis:polyOrder()

   if self.isElectromagnetic then
      -- solve for initial Apar
      local apar = self.potentials[1].apar
      self.currentDens:clear(0.0)
      for nm, s in pairs(species) do
         self.currentDens:accumulate(s:getCharge(), s:getMomDensity())
      end
      self.aparSlvr:advance(0.0, {self.currentDens}, {apar})

      -- decrease effective polynomial order in z of apar by setting the highest order z coefficients to 0
--      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
--         local localRange = apar:localRange()
--         local indexer = apar:genIndexer()
--         local ptr = apar:get(1)
--
--         -- loop over all cells
--         for idx in localRange:rowMajorIter() do
--            self.grid:setIndex(idx)
--            
--            apar:fill(indexer(idx), ptr)
--            if self.ndim == 1 then
--               ptr:data()[polyOrder] = 0.0
--            else -- ndim == 3
--               if polyOrder == 1 then
--                  ptr:data()[3] = 0.0
--                  ptr:data()[5] = 0.0
--                  ptr:data()[6] = 0.0
--                  ptr:data()[7] = 0.0
--               elseif polyOrder == 2 then
--                  ptr:data()[9] = 0.0
--                  ptr:data()[13] = 0.0
--                  ptr:data()[14] = 0.0
--                  ptr:data()[15] = 0.0
--                  ptr:data()[16] = 0.0
--                  ptr:data()[17] = 0.0
--                  ptr:data()[18] = 0.0
--                  ptr:data()[19] = 0.0
--               end
--            end
--         end
--      end

      -- clear dApar/dt ... will be solved for before being used
      self.potentials[1].dApardt:clear(0.0)
   end

   -- apply BCs and update ghosts
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
     -- include jacobian factor in metric coefficients if using linearized polarization density
     if self.linearizedPolarization then
        gxx = funcField.geo.gxxJ
        gxy = funcField.geo.gxyJ
        gyy = funcField.geo.gyyJ
     else  -- if not, jacobian already included in polarization density
        gxx = funcField.geo.gxx
        gxy = funcField.geo.gxy
        gyy = funcField.geo.gyy
     end
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

      self.laplacianWeight:combine(laplacianConstant, self.unitWeight) -- no jacobian here
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
       gxx = funcField.geo.gxxJ,
       gxy = funcField.geo.gxyJ,
       gyy = funcField.geo.gyyJ,
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
       gxx = funcField.geo.gxxJ,
       gxy = funcField.geo.gxyJ,
       gyy = funcField.geo.gyyJ,
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

     -- separate solver for additional step for p=1
     if self.basis:polyOrder() == 1 then
        self.dApardtSlvr2 = Updater.FemPoisson {
          onGrid = self.grid,
          basis = self.basis,
          bcLeft = self.aparBcLeft,
          bcRight = self.aparBcRight,
          bcBottom = self.aparBcBottom,
          bcTop = self.aparBcTop,
          periodicDirs = self.periodicDirs,
          zContinuous = not self.discontinuousApar,
          gxx = funcField.geo.gxxJ,
          gxy = funcField.geo.gxyJ,
          gyy = funcField.geo.gyyJ,
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
        if laplacianConstant ~= 0 then self.dApardtSlvr2:setLaplacianWeight(self.laplacianWeight) end
        if modifierConstant ~= 0 then self.dApardtSlvr2:setModifierWeight(self.modifierWeight) end
     end
   end

   -- need to set this flag so that field calculated self-consistently at end of full RK timestep
   self.isElliptic = true

   if self.isElectromagnetic and self.basis:polyOrder() == 1 then 
      self.nstep = 3 
   elseif self.isElectromagnetic then
      self.nstep = 2
   end

   -- function to construct a BC updater
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
	 onGrid = self.grid,
	 boundaryConditions = bcList,
	 dir = dir,
	 edge = edge,
      }
   end

   -- various functions to apply BCs of different types
   local function bcOpen(dir, tm, idxIn, fIn, fOut)
      -- Requires skinLoop = "pointwise".
      self.basis:flipSign(dir, fIn, fOut)
   end

   -- functions to make life easier while reading in BCs to apply
   self.boundaryConditions = { } -- list of Bcs to apply
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
   
   -- for non-periodic dirs, use BC_OPEN to make sure values on edge of ghost cells match values on edge of skin cells, 
   -- so that field is continuous across skin-ghost boundary
   for dir = 1, self.ndim do
      if not contains(self.periodicDirs, dir) then 
         handleBc(dir, {EM_BC_OPEN, EM_BC_OPEN}) 
      end
   end
end

function GkField:createDiagnostics()
   -- create Adios object for field I/O
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.potentials[1].phi:elemType(),
      method = self.ioMethod,
      writeGhost = self.writeGhost,
      metaData = {
	 polyOrder = self.basis:polyOrder(),
	 basisType = self.basis:id()
      },
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
      if self.linearizedPolarization then
         local esEnergyFac = .5*self.polarizationWeight
         if self.ndim == 1 then esEnergyFac = esEnergyFac*self.kperp2 end
         if self.energyCalc then self.energyCalc:advance(tm, { self.potentials[1].phi, esEnergyFac }, { self.esEnergy }) end
      else
         -- something
      end
      if self.isElectromagnetic then 
        local emEnergyFac = .5/self.mu0
        if self.ndim == 1 then emEnergyFac = emEnergyFac*self.kperp2 end
        if self.energyCalc then self.energyCalc:advance(tm, { self.potentials[1].apar, emEnergyFac}, { self.emEnergy }) end
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
      -- make sure phi is continuous across skin-ghost boundary
      for _, bc in ipairs(self.boundaryConditions) do
         bc:advance(tCurr, {}, {potCurr.phi})
      end
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
--      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
--         local localRange = dApardt:localRange()
--         local indexer = dApardt:genIndexer()
--         local ptr = dApardt:get(1)
--
--         -- loop over all cells
--         for idx in localRange:rowMajorIter() do
--            self.grid:setIndex(idx)
--            
--            dApardt:fill(indexer(idx), ptr)
--            if self.ndim == 1 then
--               ptr:data()[polyOrder] = 0.0
--            else -- ndim == 3
--               if polyOrder == 1 then
--                  ptr:data()[3] = 0.0
--                  ptr:data()[5] = 0.0
--                  ptr:data()[6] = 0.0
--                  ptr:data()[7] = 0.0
--               elseif polyOrder == 2 then
--                  ptr:data()[9] = 0.0
--                  ptr:data()[13] = 0.0
--                  ptr:data()[14] = 0.0
--                  ptr:data()[15] = 0.0
--                  ptr:data()[16] = 0.0
--                  ptr:data()[17] = 0.0
--                  ptr:data()[18] = 0.0
--                  ptr:data()[19] = 0.0
--               end
--            end
--         end
--      end

      -- apply BCs
      local tmStart = Time.clock()
      -- make sure dApardt is continuous across skin-ghost boundary
      for _, bc in ipairs(self.boundaryConditions) do
         bc:advance(tCurr, {}, {potCurr.dApardt})
      end
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
      self.dApardtSlvr2:setModifierWeight(self.modifierWeight)
      -- dApar/dt solve
      local dApardt = potCurr.dApardt
      self.dApardtSlvr2:advance(tCurr, {self.currentDens}, {dApardt}) 
      
      -- decrease effective polynomial order in z of dApar/dt by setting the highest order z coefficients to 0
      -- this ensures that dApar/dt is in the same space as dPhi/dz
--      if self.ndim == 1 or self.ndim == 3 then -- only have z direction in 1d or 3d (2d is assumed to be x,y)
--         local localRange = dApardt:localRange()
--         local indexer = dApardt:genIndexer()
--         local ptr = dApardt:get(1)
--
--         -- loop over all cells
--         for idx in localRange:rowMajorIter() do
--            self.grid:setIndex(idx)
--            
--            dApardt:fill(indexer(idx), ptr)
--            if self.ndim == 1 then
--               ptr:data()[polyOrder] = 0.0
--            else -- ndim == 3
--               if polyOrder == 1 then
--                  ptr:data()[3] = 0.0
--                  ptr:data()[5] = 0.0
--                  ptr:data()[6] = 0.0
--                  ptr:data()[7] = 0.0
--               elseif polyOrder == 2 then
--                  ptr:data()[9] = 0.0
--                  ptr:data()[13] = 0.0
--                  ptr:data()[14] = 0.0
--                  ptr:data()[15] = 0.0
--                  ptr:data()[16] = 0.0
--                  ptr:data()[17] = 0.0
--                  ptr:data()[18] = 0.0
--                  ptr:data()[19] = 0.0
--               end
--            end
--         end
--      end

      -- apply BCs
      local tmStart = Time.clock()
      -- make sure dApardt is continuous across skin-ghost boundary
      for _, bc in ipairs(self.boundaryConditions) do
         bc:advance(tCurr, {}, {potCurr.dApardt})
      end
      dApardt:sync(true)
      self.bcTime = self.bcTime + (Time.clock()-tmStart)

      -- Apar RHS is just dApar/dt
      potRhs.apar:copy(dApardt)
   end
end

function GkField:applyBcIdx(tCurr, idx)
   -- don't do anything here. BCs handled in advance steps.
end

-- NOTE: global boundary conditions handled by solver. this just updates interproc ghosts.
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
      syncPeriodicDirs = false,
   }

   -- bmagInv = 1/B
   self.geo.bmagInv = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- gradpar = J B / sqrt(g_zz)
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

   -- total jacobian, including phase space jacobian
   -- jacobTot = J B 
   self.geo.jacobTot = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- inverse of total jacobian, including phase space jacobian
   -- jacobTotInv = 1 / ( J B )
   self.geo.jacobTotInv = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }

   -- functions for magnetic drifts 
   -- geoX = g_xz / ( sqrt(g_zz) )
   self.geo.geoX = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- geoY = g_yz / ( sqrt(g_zz) )
   self.geo.geoY = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   -- geoZ = sqrt(g_zz)
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

   -- functions for laplacian, including jacobian factor
   self.geo.gxxJ = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   self.geo.gxyJ = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {1, 1},
      syncPeriodicDirs = false
   }
   self.geo.gyyJ = DataStruct.Field {
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
      metaData = {
	 polyOrder = self.basis:polyOrder(),
	 basisType = self.basis:id()
      },
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
      self.gxxJ_Func = function (t, xn)
         return self.jacobGeoFunc(t,xn)
      end
      self.gxyJ_Func = function (t, xn)
         return 0
      end
      self.gyyJ_Func = function (t, xn)
         return self.jacobGeoFunc(t,xn)
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
      self.gxxJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[1]*self.jacobGeoFunc(t,xn)
      end
      self.gxyJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[2]*self.jacobGeoFunc(t,xn)
      end
      self.gyyJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[3]*self.jacobGeoFunc(t,xn)
      end
   elseif self.ndim == 3 then
      self.jacobGeoFunc = function (t, xn)
         return self.grid:calcJacobian(xn)
      end
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
      self.gxxJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[1]*self.jacobGeoFunc(t,xn)
      end
      self.gxyJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[2]*self.jacobGeoFunc(t,xn)
      end
      self.gyyJ_Func = function (t, xn)
         local g = {}
         self.grid:calcContraMetric(xn, g)
         return g[4]*self.jacobGeoFunc(t,xn)
      end
   end

   self.bmagInvFunc = function (t, xn)
      return 1.0/self.bmagFunc(t,xn)
   end
   self.gradparFunc = function (t, xn)
      return self.jacobGeoFunc(t,xn)*self.bmagFunc(t,xn)/math.sqrt(self.g_zzFunc(t,xn))
   end
   self.geoXFunc = function (t, xn)
      return self.g_xzFunc(t,xn)/math.sqrt(self.g_zzFunc(t,xn))
   end
   self.geoYFunc = function (t, xn)
      return self.g_yzFunc(t,xn)/math.sqrt(self.g_zzFunc(t,xn))
   end
   self.geoZFunc = function (t, xn)
      return math.sqrt(self.g_zzFunc(t,xn))
   end

   -- inverse of jacobGeo, for removing jacobGeo factor from output quantities, e.g. density
   self.jacobGeoInvFunc = function (t, xn)
      return 1.0/self.jacobGeoFunc(t, xn)
   end

   -- total jacobian, including geo jacobian and phase jacobian
   self.jacobTotFunc = function (t, xn)
      return self.jacobGeoFunc(t, xn)*self.bmagFunc(t, xn)
   end

   -- inverse of total jacobian, including geo jacobian and phase jacobian
   self.jacobTotInvFunc = function (t, xn)
      return 1.0/(self.jacobGeoFunc(t, xn)*self.bmagFunc(t, xn))
   end

   -- projection updaters
   self.setBmag = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.bmagFunc,
      projectOnGhosts = true,
   }
   self.setBmagInv = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.bmagInvFunc
   }
   self.setGradpar = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gradparFunc
   }
   self.setJacobGeo = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobGeoFunc
   }
   self.setJacobGeoInv = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobGeoInvFunc
   }
   self.setJacobTot = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobTotFunc
   }
   self.setJacobTotInv = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.jacobTotInvFunc
   }
   self.setGeoX = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoXFunc
   }
   self.setGeoY = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoYFunc
   }
   self.setGeoZ = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.geoZFunc
   }
   self.setGxx = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxx_Func
   }
   self.setGxy = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxy_Func
   }
   self.setGyy = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gyy_Func
   }
   self.setGxxJ = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxxJ_Func
   }
   self.setGxyJ = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gxyJ_Func
   }
   self.setGyyJ = Updater.EvalOnNodes {
      onGrid = self.grid,
      basis = self.basis,
      projectOnGhosts = true,
      evaluate = self.gyyJ_Func
   }
   if self.phiWallFunc then 
      self.setPhiWall = Updater.EvalOnNodes {
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
   if self.ndim == 3 then self.bmagVars = {1,3} end

   self.smoother = Updater.FemPoisson {
     onGrid = self.grid,
     basis = self.basis,
     periodicDirs = self.periodicDirs,
     smooth = true,
   }

   self.weakDivision = Updater.CartFieldBinOp {
      onGrid = self.grid,
      weakBasis = self.basis,
      operation = "Divide",
      onGhosts = true,
   }

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
end

function GkGeometry:createDiagnostics()
end

function GkGeometry:initField()
   self.setBmag:advance(0.0, {}, {self.geo.bmag})
   self.setBmagInv:advance(0.0, {}, {self.geo.bmagInv})
   self.setGradpar:advance(0.0, {}, {self.geo.gradpar})
   self.setJacobGeo:advance(0.0, {}, {self.geo.jacobGeo})
   self.setJacobGeoInv:advance(0.0, {}, {self.geo.jacobGeoInv})
   self.setJacobTot:advance(0.0, {}, {self.geo.jacobTot})
   self.setJacobTotInv:advance(0.0, {}, {self.geo.jacobTotInv})
   self.setGxx:advance(0.0, {}, {self.geo.gxx})
   self.setGxy:advance(0.0, {}, {self.geo.gxy})
   self.setGyy:advance(0.0, {}, {self.geo.gyy})
   self.setGxxJ:advance(0.0, {}, {self.geo.gxxJ})
   self.setGxyJ:advance(0.0, {}, {self.geo.gxyJ})
   self.setGyyJ:advance(0.0, {}, {self.geo.gyyJ})
   self.setGeoX:advance(0.0, {}, {self.geo.geoX})
   self.setGeoY:advance(0.0, {}, {self.geo.geoY})
   self.setGeoZ:advance(0.0, {}, {self.geo.geoZ})
   if self.setPhiWall then self.setPhiWall:advance(0.0, {}, {self.geo.phiWall})
   else self.geo.phiWall:clear(0.0) end

   -- smooth all geo quantities
   -- note bmag is only one that is required to be smooth, but
   -- for others would need to pass R and L values in surface kernels if not smooth
--   self.smoother:advance(0.0, {self.geo.bmag}, {self.geo.bmag})
--   self.smoother:advance(0.0, {self.geo.jacobTotInv}, {self.geo.jacobTotInv})
--   --self.smoother:advance(0.0, {self.geo.bmagInv}, {self.geo.bmagInv})
--   self.smoother:advance(0.0, {self.geo.gradpar}, {self.geo.gradpar})
--   self.smoother:advance(0.0, {self.geo.jacobGeo}, {self.geo.jacobGeo})
--   --self.smoother:advance(0.0, {self.geo.jacobGeoInv}, {self.geo.jacobGeoInv})
--   --self.smoother:advance(0.0, {self.geo.jacobTotInv}, {self.geo.jacobTotInv})
--   self.smoother:advance(0.0, {self.geo.gxx}, {self.geo.gxx})
--   self.smoother:advance(0.0, {self.geo.gxy}, {self.geo.gxy})
--   self.smoother:advance(0.0, {self.geo.gyy}, {self.geo.gyy})
--   self.smoother:advance(0.0, {self.geo.geoX}, {self.geo.geoX})
--   self.smoother:advance(0.0, {self.geo.geoY}, {self.geo.geoY})
--   self.smoother:advance(0.0, {self.geo.geoZ}, {self.geo.geoZ})

   --self.weakDivision:advance(0.0, {self.geo.jacobTot, self.unitWeight}, {self.geo.jacobTotInv})
   --self.weakDivision:advance(0.0, {self.geo.jacobGeo, self.unitWeight}, {self.geo.jacobGeoInv})

   -- sync ghost cells. these calls do not enforce periodicity because
   -- these fields initialized with syncPeriodicDirs = false
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.gradpar:sync(false)
   self.geo.gxx:sync(false)
   self.geo.gxy:sync(false)
   self.geo.gyy:sync(false)
   self.geo.gxxJ:sync(false)
   self.geo.gxyJ:sync(false)
   self.geo.gyyJ:sync(false)
   self.geo.jacobGeo:sync(false)
   self.geo.jacobGeoInv:sync(false)
   self.geo.jacobTotInv:sync(false)
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
      self.fieldIo:write(self.geo.gradpar, string.format("gradpar_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoX, string.format("geoX_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoY, string.format("geoY_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.geoZ, string.format("geoZ_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gxx, string.format("gxx_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gxy, string.format("gxy_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.gyy, string.format("gyy_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.jacobGeo, string.format("jacobGeo_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.fieldIo:write(self.geo.jacobTotInv, string.format("jacobTotInv_%d.bp", self.ioFrame), tm, self.ioFrame)
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
