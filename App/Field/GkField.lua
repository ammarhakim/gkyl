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
local Range            = require "Lib.Range"
local DecompRegionCalc = require "Lib.CartDecomp"
local diff             = require("sci.diff")
local FieldBase        = require "App.Field.FieldBase"
local Grid             = require "Grid"
local LinearTrigger    = require "LinearTrigger"
local Logger           = require "Lib.Logger"
local lume             = require "Lib.lume"
local math             = require("sci.math").generic
local Proto            = require "Lib.Proto"
local Species          = require "App.Species"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local xsys             = require "xsys"

local GkField = Proto(FieldBase.FieldBase)

local checkBCs = function(dimIn, isDirPer, bcLoIn, bcUpIn, bcLoOut, bcUpOut)
   local periodicDomain = true
   for d=1,dimIn do periodicDomain = periodicDomain and isDirPer[d] end
   if bcLoIn==nil and bcUpIn==nil then
      if periodicDomain then
         for d=1,dimIn do bcLoOut[d], bcUpOut[d] = {T="P"}, {T="P"} end
      else
         assert(dimIn==1, "App.Field.GkField: must specify 'bcLower' and 'bcUpper' if dimensions > 1.")
      end
   else
      assert(#bcLoIn==#bcUpIn, "App.Field.GkField: number of entries in bcLower and bcUpper must be equal.")
      assert(dimIn==1 or (dimIn>1 and #bcLoIn>=2), "App.Field.GkField: number of entries in bcLower/bcUpper must >= 2.")
      for d=1,#bcLoIn do bcLoOut[d], bcUpOut[d] = bcLoIn[d], bcUpIn[d] end
   end
   for d=1,math.min(dimIn,2) do
      if isDirPer[d] then
         assert(bcLoOut[d].T=="P" and bcUpOut[d].T=="P",
                string.format("App.Field.GkField: direction %d is periodic. Must use {T='P'} in bcLower/bcUpper.",d))
     
      end
   end
end


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
   
   self.evolve = xsys.pickBool(tbl.evolve, true) -- By default evolve field.

   self.isElectromagnetic = xsys.pickBool(tbl.isElectromagnetic, false) -- Electrostatic by default.

   -- Create triggers to write fields.
   local nFrame = tbl.nFrame or appTbl.nFrame
   self.ioTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   self.ioFrame = 0 -- Frame number for IO.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   
   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   self.externalPhi = tbl.externalPhi
   if (self.externalPhi and self.evolve) then print("GkField: will not solve Poisson problem and use external phi instead.") end

   -- This allows us to apply a multiplicative time dependence to externalPhi.
   if tbl.externalPhiTimeDependence then
      self.externalPhiTimeDependence = tbl.externalPhiTimeDependence
   else
      self.externalPhiTimeDependence = function(t) return 1.0 end
   end

   if not self.externalPhi then
      -- Get boundary conditions.
      local ndim = #appTbl.lower
      if appTbl.periodicDirs then self.periodicDirs = appTbl.periodicDirs else self.periodicDirs = {} end
      local isDirPeriodic = {}
      for d = 1, ndim do isDirPeriodic[d] = lume.find(self.periodicDirs,d) ~= nil end
      self.bcLowerPhi, self.bcUpperPhi   = {}, {}
      self.bcLowerApar, self.bcUpperApar = {}, {}
   
      -- Allow unspecified BCs if domain is periodic. Or if domain is not periodic
      -- allow unspecified z-BCs, but do not allow unspecified xy-BCs for ndim>1.
      assert((tbl.bcLowerPhi and tbl.bcUpperPhi) or (tbl.bcLowerPhi==nil and tbl.bcUpperPhi==nil),
             "App.Field.GkField: must specify both 'bcLowerPhi' and 'bcUpperPhi' or none.")
      if #self.bcLowerPhi==0 and #self.bcUpperPhi==0 then  -- Needed to not override the backward compatible part above.
         checkBCs(ndim, isDirPeriodic, tbl.bcLowerPhi, tbl.bcUpperPhi, self.bcLowerPhi, self.bcUpperPhi)
      end
      if self.isElectromagnetic then
         assert((tbl.bcLowerApar and tbl.bcUpperApar) or (tbl.bcLowerApar==nil and tbl.bcUpperApar==nil),
                "App.Field.GkField: must specify both 'bcLowerApar' and 'bcUpperApar' or none.")
         if #self.bcLowerApar==0 and #self.bcUpperApar==0 then  -- Needed to not override the backward compatible part above.
            checkBCs(ndim, isDirPeriodic, tbl.bcLowerApar, tbl.bcUpperApar, self.bcLowerApar, self.bcUpperApar)
         end
      end

      self.adiabatic = false
      self.discontinuousPhi  = xsys.pickBool(tbl.discontinuousPhi, false)
      self.discontinuousApar = xsys.pickBool(tbl.discontinuousApar, true)

      -- For ndim=1 only.
      self.kperpSq = tbl.kperpSq or tbl.kperp2   -- kperp2 for backwards compatibility.
      assert((ndim==1 and self.kperpSq) or ndim>1, "GkField: must specify kperpSq for ndim=1")

      -- Determine whether to use linearized polarization term in poisson equation,
      -- which uses background density in polarization weight.
      -- If not, uses full time-dependent density in polarization weight.
      self.linearizedPolarization = xsys.pickBool(tbl.linearizedPolarization, true)
      self.uniformPolarization    = xsys.pickBool(tbl.uniformPolarization, self.linearizedPolarization)
      assert(not (self.linearizedPolarization==false and self.uniformPolarization==true), "GkField: cannot have nonlinearized and uniform polarization.")

      if self.isElectromagnetic then self.mu0 = tbl.mu0 or Constants.MU0 end

   end

   -- Create trigger for how frequently to compute field energy.
   -- Do not compute the integrated diagnostics less frequently than we output data.
   if appTbl.calcIntQuantEvery then
      self.calcIntFieldEnergyTrigger = LinearTrigger(0, appTbl.tEnd,  math.max(nFrame,math.floor(1/appTbl.calcIntQuantEvery)))
   else
      self.calcIntFieldEnergyTrigger = function(t) return true end
   end

   self.timers = {advance = 0.,   bc = 0.}
end

-- Methods for EM field object.
function GkField:hasEB() return true, self.isElectromagnetic end
function GkField:setGrid(grid)
   self.grid = grid
   self.ndim = self.grid:ndim()

   local keepDims = {};  for i = 1, self.ndim do keepDims[i] = i end
   local gridInfo = grid:childGrid(keepDims)
   local GridConstructor = grid.mapc2p and Grid.MappedCart or Grid.RectCart

   -- Create global in z (local in xy) grid.
   local xyComm = self.grid:getMessenger():getSubComms_host()['xy']
   local xyCuts = {};  xyCuts[self.ndim] = 1
   for d = self.ndim-1, 1, -1 do xyCuts[d] = self.grid:cuts(d) end
   local xyDecomp = DecompRegionCalc.CartProd {  -- Decomp in x-y but not z.
      cuts = xyCuts,  comm = xyComm,
   }
   self.gridGlobalZ = GridConstructor {
      lower = gridInfo.lower,  world    = gridInfo.world,
      upper = gridInfo.upper,  mappings = gridInfo.coordinateMap,
      cells = gridInfo.cells,  mapc2p   = gridInfo.mapc2p,
      periodicDirs  = gridInfo.periodicDirs,
      messenger     = gridInfo.messenger,
      decomposition = xyDecomp,
   }

   -- Create global in xy (local in z) grid.
   local zComm = self.grid:getMessenger():getSubComms_host()['z']
   local zCuts = {};  zCuts[self.ndim] = self.grid:cuts(self.ndim)
   for d = self.ndim-1, 1, -1 do zCuts[d] = 1 end
   local zDecomp = DecompRegionCalc.CartProd {  -- Decomp in x-y but not z.
      cuts = zCuts,  comm = zComm,
   }
   self.gridGlobalXY = GridConstructor {
      lower = gridInfo.lower,  world    = gridInfo.world,
      upper = gridInfo.upper,  mappings = gridInfo.coordinateMap,
      cells = gridInfo.cells,  mapc2p   = gridInfo.mapc2p,
      periodicDirs  = gridInfo.periodicDirs,
      messenger     = gridInfo.messenger,
      decomposition = zDecomp,
   }
end

local function createField(grid, basis, ghostCells, vComp, periodicSync, useDevice)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid           = grid,
      numComponents    = basis:numBasis()*vComp,
      ghost            = ghostCells,
      metaData         = {polyOrder = basis:polyOrder(),
                          basisType = basis:id()},
      syncPeriodicDirs = periodicSync,
      useDevice        = useDevice,
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
      self.potentials[i].apar    = createField(self.grid,self.basis,{1,1})
      self.potentials[i].dApardt = createField(self.grid,self.basis,{1,1})
      self.potentials[i].apar:clear(0.0)
      self.potentials[i].dApardt:clear(0.0)
   end

   self.dApardtProv = createField(self.grid,self.basis,{1,1})

   -- Create fields for total charge density and other things needed for Poisson solve.
   self.chargeDens              = createField(self.grid,self.basis,{1,1})
   self.chargeDensLocal         = createField(self.grid,self.basis,{1,1})
   self.polarizationWeight      = createField(self.grid,self.basis,{1,1})
   self.polarizationWeightLocal = createField(self.grid,self.basis,{1,1})
   self.modWeightPoisson        = createField(self.grid,self.basis,{1,1})

   if self.isElectromagnetic then
      -- Create fields for total current density and other things needed for Ampere solve.
      self.currentDens           = createField(self.grid,self.basis,{1,1})
      self.currentDensLocal      = createField(self.grid,self.basis,{1,1})
      self.dApardtSlvr_lapModFac = createField(self.grid,self.basis,{1,1})
      self.lapWeightAmpere       = createField(self.grid,self.basis,{1,1})
      self.modWeightAmpere       = createField(self.grid,self.basis,{1,1})
      self.modWeightAmpereLocal  = createField(self.grid,self.basis,{1,1})
   end

   -- Extra fields needed for the distributed parallel field solve.
   self.localBuffer    = createField(self.grid,self.basis,{0,0})
   self.globalBufferZ  = createField(self.gridGlobalZ,self.basis,{0,0})
   self.globalSolZ     = createField(self.gridGlobalZ,self.basis,{1,1})
   self.globalBufferXY = createField(self.gridGlobalXY,self.basis,{0,0})
   self.globalSolXY    = createField(self.gridGlobalXY,self.basis,{1,1})
   
   -- For diagnostics:
   self.phiSq = createField(self.grid,self.basis,{1,1})
   self.esEnergyFac = createField(self.grid,self.basis,{1,1})

   -- For storing integrated energies.
   self.intPhiSq          = DataStruct.DynVector { numComponents = 1, }
   self.esEnergyAdiabatic = DataStruct.DynVector { numComponents = 1, }
   self.gradPerpPhiSq     = DataStruct.DynVector { numComponents = 1, }
   self.aparSq            = DataStruct.DynVector { numComponents = 1, }
   self.esEnergy          = DataStruct.DynVector { numComponents = 1, }
   self.emEnergy          = DataStruct.DynVector { numComponents = 1, }
end

-- Solve for initial fields self-consistently 
-- from initial distribution function.
function GkField:initField(population)
   if not self.externalPhi then
      -- Solve for initial phi. Temporarily re-set the calcPhi function in the advance method.
      local trueCalcPhi = self.calcPhi
      self.calcPhi = (not self.evolve) and self.calcPhiInitial or self.calcPhi

      self:advance(0.0, population, 1, 1)

      if (not self.evolve) then self.calcPhi = trueCalcPhi end
   end

   if self.isElectromagnetic then
      -- Solve for initial Apar.
      local apar = self.potentials[1].apar
      self.currentDensLocal:clear(0.0)
      for _, s in population.iterLocal() do
         self.currentDensLocal:accumulate(s:getCharge(), s:getMomDensity())
      end
      -- Reduce currents across species communicator.
      self.currentDens:clear(0.0)
      population:AllreduceByCell(self.currentDensLocal, self.currentDens, 'sum')

      self.aparSlvr:advance(0.0, {self.currentDens}, {apar})

      -- Clear dApar/dt ... will be solved for before being used.
      self.potentials[1].dApardt:clear(0.0)
   end

   -- Apply BCs and update ghosts.
   self:applyBc(0, self.potentials[1])

   if self.ndim > 1 and self.ioFrame == 0 and (not self.externalPhi) then 
      self.fieldIo:write(self.modWeightPoisson,  "modifierWeightPoisson_0.bp", tm, self.ioFrame, false)
   end
end

function GkField:rkStepperFields() return self.potentials end

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

function GkField:createSolver(population, externalField)
   -- Need to set this flag so that field calculated self-consistently at end of full RK timestep.
   self.isElliptic = true

   if self.externalPhi then
      local evalOnNodes = Updater.EvalOnNodes {
         onGrid = self.grid,   evaluate = self.externalPhi,
         basis  = self.basis,  onGhosts = true
      }
      self.externalPhiFld = createField(self.grid,self.basis,{1,1})
      evalOnNodes:advance(0.0, {}, {self.externalPhiFld})
      for i = 1, self.nRkDup do
         self.potentials[i].phi:combine(self.externalPhiTimeDependence(0.0),self.externalPhiFld)
      end

      if self.evolve then
	 self.calcPhi = function(tCurr, populationIn, potentialsIn, potentialsOut)
            potentialsIn.phi:combine(self.externalPhiTimeDependence(tCurr), self.externalPhiFld)
	 end
      else
	 self.calcPhi = function(tCurr, populationIn, potentialsIn, potentialsOut) end
      end
   else
      local perpDim     = 2
      local parDir      = self.ndim == 2 and -1 or self.ndim
      local decompRange = self.grid:decomposedRange()
      local cutsRange   = decompRange:cutsRange()  -- Range of MPI processes sharing a field solve.
      local subdomIdx   = {};  decompRange:cutsInvIndexer()(self.grid:subGridId(), subdomIdx)
      -- Range of MPI processes participating in z-field solve.
      local removeDir, locInDir = {}, {}
      for d=1,perpDim do removeDir[d], locInDir[d] = 1, subdomIdx[d] end
      removeDir[parDir], locInDir[parDir] = 0, 0
      self.zCutsRange = cutsRange:deflate(removeDir, locInDir)
      -- Ranges in a z-global field into which we will copy a buffer obtained form other processes.
      self.zDestRange, self.zBufferOffset = {}, {}
      for zIdx in self.zCutsRange:colMajorIter() do
         local idx = {};  for d=1,perpDim do idx[d] = subdomIdx[d] end
         idx[parDir] = zIdx[1]

         local linIdx    = decompRange:cutsIndexer()(idx)
         local destRange = decompRange:subDomain(linIdx)
         local lv, uv    = destRange:lowerAsVec(), destRange:upperAsVec()
         self.zDestRange[zIdx[1]]    = self.globalSolZ:localExtRange():subRange(lv,uv)
         self.zBufferOffset[zIdx[1]] = (zIdx[1]-1)*self.localBuffer:size()
      end
      -- Range in a global field we'll copy the local solution from.
      local lv, uv = self.potentials[1].phi:localRange():lowerAsVec(), self.potentials[1].phi:localRange():upperAsVec()
      self.zSrcRange = self.globalSolZ:localExtRange():subRange(lv,uv)

      if self.ndim > 1 then
         -- Range of MPI processes participating in xy-field solve.
         local removeDir, locInDir = {}, {}
         for d=1,perpDim do removeDir[d], locInDir[d] = 0, 0 end
         removeDir[parDir], locInDir[parDir] = 1, subdomIdx[parDir]
         self.xyCutsRange = cutsRange:deflate(removeDir, locInDir)
         -- Ranges in a xy-global field into which we will copy a buffer obtained form other processes.
         self.xyDestRange, self.xyBufferOffset = {}, {}
         self.xyCutsRangeIdxr = self.xyCutsRange:indexer(Range.colMajor)
         for xyIdx in self.xyCutsRange:colMajorIter() do
            local linIdxPerp = self.xyCutsRangeIdxr(xyIdx[1],xyIdx[2])

            local idx = {};  for d=1,perpDim do idx[d] = xyIdx[d] end
            idx[parDir] = subdomIdx[parDir]

            local linIdx    = decompRange:cutsIndexer()(idx)
            local destRange = decompRange:subDomain(linIdx)
            local lv, uv    = destRange:lowerAsVec(), destRange:upperAsVec()
            self.xyDestRange[linIdxPerp]    = self.globalSolXY:localExtRange():subRange(lv,uv)
            self.xyBufferOffset[linIdxPerp] = (linIdxPerp-1)*self.localBuffer:size()
         end
         -- Range in a global field we'll copy the local solution from.
         local lv, uv = self.potentials[1].phi:localRange():lowerAsVec(), self.potentials[1].phi:localRange():upperAsVec()
         self.xySrcRange = self.globalSolXY:localExtRange():subRange(lv,uv)
      end

      -- Get adiabatic species info.
      for _, s in population.iterGlobal() do
         if Species.AdiabaticSpecies.is(s) then
            self.adiabatic, self.adiabSpec = true, s
         end
      end
      assert((self.adiabatic and self.isElectromagnetic) == false, "GkField: cannot use adiabatic response for electromagnetic case")

      self.unitField = createField(self.grid,self.basis,{1,1})
      self.unitField:shiftc(1.*(math.sqrt(2.)^(self.ndim)), 0)

      -- Metric coefficients in the Poisson and Ampere equations for phi and Apar.
      local gxxPoisson, gxyPoisson, gyyPoisson
      local gxxAmpere, gxyAmpere, gyyAmpere
      if externalField.geo then 
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
         self.jacobGeo = externalField.geo.jacobGeo
      else
         self.jacobGeo = self.unitField
      end

      self.weakMultiply = Updater.CartFieldBinOp {
         weakBasis = self.basis,  operation = "Multiply",
         onGhosts  = true,
      }

      -- When using a linearizedPolarization term in Poisson equation,
      -- the density in the polarization density is constant in time.
      if self.linearizedPolarization then
         -- Calculate weight on polarization term: sum_s m_s * jacobGeo * n_s / B^2.
         self.polarizationWeightLocal:clear(0.)
         for _, s in population.iterLocal() do
            if Species.GkSpecies.is(s) then
               self.polarizationWeightLocal:accumulate(1., s:getPolarizationWeight())
            end
         end
         -- Reduce polarization weight across species communicator.
         self.polarizationWeight:clear(0.0)
         population:AllreduceByCell(self.polarizationWeightLocal, self.polarizationWeight, 'sum')

         self.esEnergyFac:combine(.5, self.polarizationWeight)

         -- Set up scalar multipliers on laplacian and modifier terms in Poisson equation.
         if self.ndim==1 then
            self.modWeightPoisson:combine(self.kperpSq, self.polarizationWeight)
            self.esEnergyFac:scale(self.kperpSq)  -- For diagnostics.
            -- Add contribution from the adiabatic species.
            if self.adiabatic then 
               self.modWeightPoisson:accumulate(1., self.adiabSpec:getQneutFac())
               self.esEnergyFac:accumulate(0.5, self.adiabSpec:getQneutFac())   -- For diagnostics.
            end
         else 
            self.modWeightPoisson:clear(0.)
            if self.adiabatic then   -- Add contribution from the adiabatic species.
               self.modWeightPoisson:accumulate(1., self.adiabSpec:getQneutFac())
            end
         end

         self.setPolarizationWeightFunc = function(populationTbl) end  -- Polarization weight is not time dependent.
      else
         -- When not using linearizedPolarization, weights are set each step in advance method using these functions.
         if self.ndim == 1 then
            self.setPolarizationWeightFunc = function(populationTbl)
               self.polarizationWeightLocal:clear(0.)
               for _, s in populationTbl.iterLocal() do
                  self.polarizationWeightLocal:accumulate(self.kperpSq, s:getPolarizationMassDensity())
               end
               -- Reduce polarization weight across species communicator.
               self.polarizationWeight:clear(0.0)
               populationTbl:AllreduceByCell(self.polarizationWeightLocal, self.polarizationWeight, 'sum')
            end
         else
            self.setPolarizationWeightFunc = function(populationTbl)
               self.polarizationWeightLocal:clear(0.)
               for _, s in populationTbl.iterLocal() do
                  self.polarizationWeightLocal:accumulate(-1.0, s:getPolarizationMassDensity())
               end
               -- Reduce polarization weight across species communicator.
               self.polarizationWeight:clear(0.0)
               populationTbl:AllreduceByCell(self.polarizationWeightLocal, self.polarizationWeight, 'sum')
            end
         end
      end

      -- MF 2022/08/23: Construction of modWeightPoisson and epsPoisson took place
      -- on the device because they were assigned with :accumulate and :combine, which use the
      -- device pointer if the CartField was created with useDevice=nil. Copy them to host.
      if self.ndim == 1 then
         -- Gather the weight in the z-field solve.
         self:AllgatherFieldZ(self.modWeightPoisson, self.globalSolZ)
         self.globalSolZ:copyDeviceToHost()  -- FemParproj expects the weight on the host.

         self.phiParSlvr = Updater.FemParproj {
            onGrid = self.gridGlobalZ,  basis   = self.basis,
            weight = self.globalSolZ,   onField = self.globalSolZ,
            periodicParallelDir = lume.any(self.periodicDirs, function(t) return t==self.ndim end),
         }
         self.phiSlvr = {
            advance = function(tcurr, inFlds, outFlds)
               local qDens  = inFlds[1]
               local phiOut = outFlds[1].phi  -- OutFlds[1] is a table of potentials.

               -- Gather the charge density onto a global field.
               self:AllgatherFieldZ(qDens, self.globalSolZ)

               -- Solve linear problem.
               self.phiParSlvr:advance(tCurr, {self.globalSolZ}, {self.globalSolZ})

               -- Copy portion of the global field belonging to this process.
               phiOut:copyRangeToRange(self.globalSolZ, phiOut:localRange(), self.zSrcRange)
            end,
	    printDevDiagnostics = function() self.phiParSlvr:printDevDiagnostics() end,
         }
      else
         local epsPoisson = createField(self.grid,self.basis,{1,1},3)
         epsPoisson:combineOffset(1., gxxPoisson, 0*self.basis:numBasis(),
                                       1., gxyPoisson, 1*self.basis:numBasis(),
                                       1., gyyPoisson, 2*self.basis:numBasis())
         self.weakMultiply:advance(0., {self.polarizationWeight, epsPoisson}, {epsPoisson})
         epsPoisson:copyDeviceToHost()
         self.modWeightPoisson:copyDeviceToHost()

         -- Gather the permittivity and the k^2 multiplier.
         local AllgatherFieldXY3 = function(fldLocal, fldGlobal)
            -- Gather a CartField with 3 times more components distributed along xy onto a global-in-xy field.
            local localBuffer    = createField(self.grid,self.basis,{0,0},3)
            local globalBufferXY = createField(self.gridGlobalXY,self.basis,{0,0},3)

            fldLocal:copyRangeToBuffer(fldLocal:localRange(), localBuffer:data())
         
            -- Gather flat buffers from other processes.
            local mess = self.grid:getMessenger()
            local comm = mess:getComms()["xy"]
            mess:Allgather(localBuffer, globalBufferXY, comm) 
         
            -- Rearrange into a global field, i.e. reconstruct the multi-D array
            -- using the collection of flat buffers we gathered from each process.
            for xyIdx in self.xyCutsRange:colMajorIter() do
               local linIdxPerp = self.xyCutsRangeIdxr(xyIdx[1],xyIdx[2])
               fldGlobal:copyRangeFromBuffer(self.xyDestRange[linIdxPerp],
                  globalBufferXY:data()+3*self.xyBufferOffset[linIdxPerp])
            end
         end
         self.epsPoissonGlobal        = createField(self.gridGlobalXY,self.basis,{1,1},3)
         local modWeightPoissonGlobal = createField(self.gridGlobalXY,self.basis,{1,1})
         AllgatherFieldXY3(epsPoisson, self.epsPoissonGlobal)
         self:AllgatherFieldXY(self.modWeightPoisson, modWeightPoissonGlobal)
         -- Permittivity and k^2 modifier were initialized on the device,
         -- but FemPoissonPerp expects it on the host.
         self.epsPoissonGlobal:copyDeviceToHost()
         modWeightPoissonGlobal:copyDeviceToHost()

         self.phiPerpSlvr = self.ndim == 2
            and Updater.FemPoisson {
               onGrid  = self.gridGlobalXY,      bcLower = self.bcLowerPhi,
               basis   = self.basis,             bcUpper = self.bcUpperPhi,
               epsilon = self.epsPoissonGlobal,  kSq     = modWeightPoissonGlobal,
            }
            or Updater.FemPoissonPerp {
               onGrid  = self.gridGlobalXY,      bcLower = self.bcLowerPhi,
               basis   = self.basis,             bcUpper = self.bcUpperPhi,
               epsilon = self.epsPoissonGlobal,  kSq     = modWeightPoissonGlobal,
            }
         self.phiPerpSolve = function(tCurr, qDensIn, phiOut)
            -- Gather the discontinuous field onto a global field.
            self:AllgatherFieldXY(qDensIn, self.globalSolXY)

            self.phiPerpSlvr:advance(tCurr, {self.globalSolXY}, {self.globalSolXY})

            -- Copy portion of the global field belonging to this process.
            phiOut:copyRangeToRange(self.globalSolXY, phiOut:localRange(), self.xySrcRange)
         end

         self.parSmooth = function(tCurr, fldIn, fldOut) fldOut:copy(fldIn) end
         if self.ndim == 3 and not self.discontinuousPhi then
            if (self.bcLowerPhi[3] and self.bcUpperPhi[3]) then
               assert(self.bcLowerPhi[3].T == "AxisymmetricLimitedTokamak", "App.Field.GkField: z BC not supported. Use a supported option or remove this z BC from the input file.")
               assert(self.bcLowerPhi[3].T == self.bcUpperPhi[3].T, "App.Field.GkField: 'bcLowerPhi[3]' and 'bcUpperPhi[3]' must be equal.")
               assert(self.bcLowerPhi[3].xLCFS == self.bcUpperPhi[3].xLCFS, "App.Field.GkField: 'bcLowerPhi[3]' and 'bcUpperPhi[3]' must be equal.")

               -- Reduce range to the core part and SOL part.
               -- Assume the split happens at a cell boundary and within the domain.
               local xLCFS = self.bcLowerPhi[3].xLCFS 
               assert(self.grid:lower(1)<xLCFS and xLCFS<self.grid:upper(1), "App.Field.GkField: 'xLCFS' coordinate must be within the x-domain.")
               local needint = (xLCFS-self.grid:lower(1))/self.grid:dx(1)
               assert(math.floor(math.abs(needint-math.floor(needint))) < 1., "App.Field.GkField: 'xLCFS' must fall on a cell boundary along x.")
               -- Determine the index of the cell that abuts xLCFS from below.
               local coordLCFS, idxLCFS = {xLCFS-1.e-7}, {-9}
               local xGridIngr = self.grid:childGrid({1})
               local xGrid = Grid.RectCart {
                  lower = xGridIngr.lower,  periodicDirs  = xGridIngr.periodicDirs,
                  upper = xGridIngr.upper,  decomposition = xGridIngr.decomposition,
                  cells = xGridIngr.cells,
               }
               xGrid:findCell(coordLCFS, idxLCFS)
               local solRange     = self.globalSolZ:localRange():shortenFromBelow(1, self.grid:numCells(1)-idxLCFS[1]+1)
               local coreRange    = self.globalSolZ:localRange():shorten(1, idxLCFS[1]+1)
               local solRangeExt  = self.globalSolZ:localExtRange():shortenFromBelow(1, self.grid:numCells(1)-idxLCFS[1]+1)
               local coreRangeExt = self.globalSolZ:localExtRange():shorten(1, idxLCFS[1]+1)
               self.solRange     = self.globalSolZ:localExtRange():subRange(    solRange:lowerAsVec(),    solRange:upperAsVec())
               self.coreRange    = self.globalSolZ:localExtRange():subRange(   coreRange:lowerAsVec(),   coreRange:upperAsVec())
               self.solRangeExt  = self.globalSolZ:localExtRange():subRange( solRangeExt:lowerAsVec(), solRangeExt:upperAsVec())
               self.coreRangeExt = self.globalSolZ:localExtRange():subRange(coreRangeExt:lowerAsVec(),coreRangeExt:upperAsVec())

               self.parSmootherCore = Updater.FemParproj {
                  onGrid  = self.gridGlobalZ,  basis               = self.basis,
                  onField = self.globalSolZ,   periodicParallelDir = true,
                  onRange = self.coreRange,    onExtRange          = self.coreRangeExt,
               }
               self.parSmootherSOL = Updater.FemParproj {
                  onGrid  = self.gridGlobalZ,  basis               = self.basis,
                  onField = self.globalSolZ,   periodicParallelDir = false,
                  onRange = self.solRange,     onExtRange          = self.solRangeExt,
               }
               self.parSmoother = {
                  advance = function(tCurr, fldIn, fldOut)
                     self.parSmootherCore:advance(tCurr, {self.globalSolZ}, {self.globalSolZ})
                     self.parSmootherSOL:advance(tCurr, {self.globalSolZ}, {self.globalSolZ})
                  end,
               }
            else
               self.parSmootherUpd = Updater.FemParproj {
                  onGrid  = self.gridGlobalZ,  basis = self.basis,
                  onField = self.globalSolZ,
                  periodicParallelDir = lume.any(self.periodicDirs, function(t) return t==self.ndim end),
               }
               self.parSmoother = {
                  advance = function(tCurr, fldIn, fldOut)
                     self.parSmootherUpd:advance(tCurr, {self.globalSolZ}, {self.globalSolZ})
                  end,
               }
            end
            self.parSmooth = function(tCurr, fldIn, fldOut)
               -- Gather the discontinuous field onto a global field.
               self:AllgatherFieldZ(fldIn, self.globalSolZ)

               self.parSmoother.advance(tCurr, {self.globalSolZ}, {self.globalSolZ})

               -- Copy portion of the global field belonging to this process.
               fldOut:copyRangeToRange(self.globalSolZ, fldOut:localRange(), self.zSrcRange)
            end
         end

         self.phiSlvr = {
            advance = function(tCurr, inFlds, outFlds)
               local qDens = inFlds[1]
               local phiAux, phiOut = outFlds[1].phiAux, outFlds[1].phi  -- OutFlds[1] is a table of potentials.
               self.parSmooth(tCurr, qDens, qDens) -- Reusing chargeDensLocal temporarily.
               self.phiPerpSolve(tCurr, qDens, phiAux)
               self.parSmooth(tCurr, phiAux, phiOut)
            end,
            printDevDiagnostics = function()
--               self.parSmoother:printDevDiagnostics()
--               self.phiPerpSlvr:printDevDiagnostics()
            end,
         }
      end

      -- Function called in the advance method.
      self.calcPhiInitial = function(tCurr, populationIn, potentialsIn, potentialsOut)
         self.chargeDensLocal:clear(0.0)
         for _, s in populationIn.iterLocal() do
            self.chargeDensLocal:accumulate(s:getCharge(), s:getNumDensity())
         end
         -- Reduce charge density across species communicator.
         self.chargeDens:clear(0.0)
         populationIn:AllreduceByCell(self.chargeDensLocal, self.chargeDens, 'sum')

         -- If not using linearized polarization term, set up laplacian weight.
         self:setPolarizationWeight(populationIn)

         -- Phi solve (elliptic, so update potentials).
         -- Energy conservation requires phi to be continuous in all directions. 
         -- The first FEM solve ensures that phi is continuous in x and y.
         -- The conserved energy is defined in terms of this intermediate result,
         -- which we denote phiAux, before the final smoothing operation in z.
         self.phiSlvr.advance(tCurr, {self.chargeDens}, {potentialsIn})
      end
      self.calcPhi = self.evolve
         and self.calcPhiInitial
	 or (self.isElectromagnetic
	     and function(tCurr, populationIn, potentialsIn, potentialsOut)
                    -- Just copy stuff over.
                    potentialsOut.apar:copy(potentialsIn.apar)
	         end
             or function(tCurr, populationIn, potentialsIn, potentialsOut) end)

   end -- end if (self.externalPhi)

   if self.isElectromagnetic then 
      local ndim = self.ndim
      local laplacianConstant, modifierConstant
      -- NOTE: aparSlvr only used to solve for initial Apar
      -- at all other times Apar is timestepped using dApar/dt
      self.aparSlvr = Updater.GkFemPoisson {
         onGrid = self.grid,   bcLower = self.bcLowerApar,
         basis  = self.basis,  bcUpper = self.bcUpperApar,
         zContinuous = not self.discontinuousApar,
         gxx = gxxAmpere,
         gxy = gxyAmpere,
         gyy = gyyAmpere,
      }
      if ndim==1 then
         laplacianConstant = 0.0
         modifierConstant  = self.kperpSq/self.mu0
      else
         laplacianConstant = -1.0/self.mu0
         modifierConstant  = 0.0
      end
      self.lapWeightAmpere:combine(laplacianConstant, self.unitField)
      self.modWeightAmpere:combine(modifierConstant, self.unitField)
      if laplacianConstant ~= 0 then self.aparSlvr:setLaplacianWeight(self.lapWeightAmpere) end
      if modifierConstant ~= 0 then self.aparSlvr:setModifierWeight(self.modWeightAmpere) end

      self.dApardtSlvr = Updater.GkFemPoisson {
         onGrid = self.grid,   bcLower = self.bcLowerApar,
         basis  = self.basis,  bcUpper = self.bcUpperApar,
         zContinuous = not self.discontinuousApar,
         gxx = gxxAmpere,
         gxy = gxyAmpere,
         gyy = gyyAmpere,
      }
      if ndim==1 then
         laplacianConstant = 0.0
         modifierConstant  = self.kperpSq/self.mu0
         self.dApardtSlvr_lapModFac:combine(modifierConstant, self.unitField) 
      else
         laplacianConstant = -1.0/self.mu0
         modifierConstant  = 1.0
         self.dApardtSlvr_lapModFac:clear(0.)
      end
      self.lapWeightAmpere:combine(laplacianConstant, self.unitField)
      self.modWeightAmpere:combine(modifierConstant, self.unitField)
      if laplacianConstant ~= 0 then self.dApardtSlvr:setLaplacianWeight(self.lapWeightAmpere) end
      if modifierConstant ~= 0 then self.dApardtSlvr:setModifierWeight(self.modWeightAmpere) end
   end

   if self.isElectromagnetic then self.nstep = 2 end

   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.potentials[1].phi:elemType(),
      writeGhost = self.writeGhost,
      metaData   = {polyOrder = self.basis:polyOrder(),
                    basisType = self.basis:id(),},
      writeRankInComm = {0, population:getComm_host(),}
   }
end

function GkField:AllgatherFieldZ(fldLocal, fldGlobal)
   -- Gather a CartField distributed along z onto a global-in-z field.
   fldLocal:copyRangeToBuffer(fldLocal:localRange(), self.localBuffer:data())

   -- Gather flat buffers from other processes.
   local mess = self.grid:getMessenger()
   local comm = mess:getComms()["z"]
   mess:Allgather(self.localBuffer, self.globalBufferZ, comm) 

   -- Rearrange into a global field, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   for zIdx in self.zCutsRange:colMajorIter() do
      fldGlobal:copyRangeFromBuffer(self.zDestRange[zIdx[1]],
         self.globalBufferZ:data()+self.zBufferOffset[zIdx[1]])
   end
end

function GkField:AllgatherFieldXY(fldLocal, fldGlobal)
   -- Gather a CartField distributed along xy onto a global-in-xy field.
   fldLocal:copyRangeToBuffer(fldLocal:localRange(), self.localBuffer:data())

   -- Gather flat buffers from other processes.
   local mess = self.grid:getMessenger()
   local comm = mess:getComms()["xy"]
   mess:Allgather(self.localBuffer, self.globalBufferXY, comm) 

   -- Rearrange into a global field, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   for xyIdx in self.xyCutsRange:colMajorIter() do
      local linIdxPerp = self.xyCutsRangeIdxr(xyIdx[1],xyIdx[2])
      fldGlobal:copyRangeFromBuffer(self.xyDestRange[linIdxPerp],
         self.globalBufferXY:data()+self.xyBufferOffset[linIdxPerp])
   end
end

function GkField:createDiagnostics()
   -- Updaters for computing integrated quantities.
   self.intCalc = Updater.CartFieldIntegrate {
      onGrid = self.grid,  basis = self.basis,
   }
   self.intSqCalc = Updater.CartFieldIntegrate {
      onGrid = self.grid,   operator = "sq",
      basis  = self.basis,
   }
   if self.ndim == 1 then
      self.energyCalc = self.intCalc
   else
      self.gradPerpSqIntCalc = Updater.CartFieldIntegrate {
         onGrid = self.grid,   operator = "gradperp_sq",
         basis  = self.basis,
      }
      self.energyCalc = Updater.CartFieldIntegrate {
         onGrid = self.grid,   operator = "eps_gradperp_sq",
         basis  = self.basis,
      }
   end

   -- The following functions streamline diagnostic computation and writing inside the time
   -- loop. This is akin to what is done in SpeciesDiagnostis, and perhaps we should use a
   -- diagnostics App for the fields too.
   self.calcPhiGridDiags = self.evolve
   and function(tm, force) end
   or function(tm, force) end
   
   self.writePhiGridDiags = self.evolve
   and function(tm, force)
          self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
       end
   or function(tm, force)
         if self.ioFrame == 0 then
            self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
         end
      end
   
   self.calcPhiIntDiags =  (self.evolve and (not self.externalPhi))
   and function(tm, force)
          -- Compute integrated quantities over domain.
          self.intSqCalc:advance(tm, { self.potentials[1].phi }, { self.intPhiSq })
          if self.energyCalc then 
             if self.ndim > 1 then
--                self.energyCalc:advance(tm, { self.potentials[1].phi, self.unitField }, { self.gradPerpPhiSq })
                self.gradPerpSqIntCalc:advance(tm, { self.potentials[1].phi }, { self.gradPerpPhiSq })
             end
             if self.linearizedPolarization then
                if self.uniformPolarization then
                   if self.ndim == 1 then 
                      self.weakMultiply:advance(0., {self.potentials[1].phi, self.potentials[1].phi},
                                                {self.phiSq})
                      self.weakMultiply:advance(0., {self.phiSq, self.esEnergyFac}, {self.phiSq})
                      self.energyCalc:advance(tm, { self.phiSq }, { self.esEnergy })
                   else
                      self.energyCalc:advance(tm, { self.potentials[1].phi, 1., self.esEnergyFac }, { self.esEnergy })
                   end
   
                   if self.adiabatic and self.ndim > 1 then
                      self.weakMultiply:advance(0., {self.potentials[1].phi, self.potentials[1].phi},
                                                {self.phiSq})
                      self.weakMultiply:advance(0., {self.phiSq, self.adiabSpec:getQneutFac()}, {self.phiSq})
                      self.intCalc:advance(tm, { self.phiSq, 0.5 }, { self.esEnergyAdiabatic })
                      local tm, energyVal = self.esEnergy:lastData()
                      local _, energyValAdiabatic = self.esEnergyAdiabatic:lastData()
                      energyVal[1] = energyVal[1] + energyValAdiabatic[1]
   --                   -- MF 2022/07/20: commenting this out for now to facilitate reg test comparison.
   --                   self.esEnergy:assignLastVal(1., energyVal)  -- Adiabatic species contribution.
                   end
                else
                   -- Something.
                end
             else
                -- Something.
             end
          end
       end
   or function(tm, force) end
   
   self.writePhiIntDiags = (self.evolve and (not self.externalPhi))
   and function(tm, force)
          self.intPhiSq:write(string.format("intPhiSq.bp"), tm, self.ioFrame)
          self.gradPerpPhiSq:write(string.format("gradPerpPhiSq.bp"), tm, self.ioFrame)
          self.esEnergy:write(string.format("esEnergy.bp"), tm, self.ioFrame)
       end
   or function(tm, force) end
   
   self.calcAparGridDiags  = function(tm, force) end
   self.writeAparGridDiags = function(tm, force) end
   self.calcAparIntDiags   = function(tm, force) end
   self.writeAparIntDiags  = function(tm, force) end

   if self.isElectromagnetic then
      self.calcAparGridDiags = self.evolve
         and function(tm, force) end
         or function(tm, force) end
      
      self.writeAparGridDiags = self.evolve
         and function(tm, force)
                self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
                self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
             end
         or function(tm, force)
               if self.ioFrame == 0 then
                  self.fieldIo:write(self.potentials[1].apar, string.format("apar_%d.bp", self.ioFrame), tm, self.ioFrame)
                  self.fieldIo:write(self.potentials[1].dApardt, string.format("dApardt_%d.bp", self.ioFrame), tm, self.ioFrame)
               end
            end
      
      local emEnergyIntCalc = self.ndim == 1 and self.intSqCalc or self.gradPerpSqIntCalc
      local emEnergyFac = .5/self.mu0
      if self.ndim == 1 then emEnergyFac = emEnergyFac*self.kperpSq end
      self.calcAparIntDiags = self.evolve
         and function(tm, force)
                -- Compute integrated quantities over domain.
                self.intSqCalc:advance(tm, { self.potentials[1].apar }, { self.aparSq })
                emEnergyIntCalc:advance(tm, { self.potentials[1].apar, emEnergyFac}, { self.emEnergy })
             end
         or function(tm, force) end
      
      self.writeAparIntDiags = self.evolve
         and function(tm, force)
                self.aparSq:write(string.format("aparSq.bp"), tm, self.ioFrame)
                self.emEnergy:write(string.format("emEnergy.bp"), tm, self.ioFrame)
             end
         or function(tm, force) end
   end

end

function GkField:write(tm, force)
   if self.calcIntFieldEnergyTrigger(tm) then
      self.calcPhiIntDiags(tm, force)
      self.calcAparIntDiags(tm, force)
   end

   if self.ioTrigger(tm) or force then
      self.writePhiGridDiags(tm, force)
      self.writeAparGridDiags(tm, force)

      self.writePhiIntDiags(tm, force)
      self.writeAparIntDiags(tm, force)
      
      self.ioFrame = self.ioFrame+1
   end
end

function GkField:writeRestart(tm)
   -- (the final "false" prevents writing of ghost cells).
   self.fieldIo:write(self.potentials[1].phi, "phi_restart.bp", tm, self.ioFrame, false)
   if self.isElectromagnetic then
      self.fieldIo:write(self.potentials[1].apar, "apar_restart.bp", tm, self.ioFrame, false)
   end
end

function GkField:readRestart()
   -- This read of restart file of phi is only to get frame
   -- numbering correct. The forward Euler recomputes the potential
   -- before updating the hyperbolic part.
   local tm, fr = self.fieldIo:read(self.potentials[1].phi, "phi_restart.bp")

   if self.isElectromagnetic then
      self.fieldIo:read(self.potentials[1].apar, "apar_restart.bp")
   end

   self:applyBc(0, self.potentials[1])

   self.ioFrame = fr 
   -- Iterate triggers.
   self.ioTrigger(tm)
end

function GkField:setPolarizationWeight(populationIn)
   self.setPolarizationWeightFunc(populationIn)
end

-- Solve for electrostatic potential phi.
function GkField:advance(tCurr, population, inIdx, outIdx)
   local tmStart = Time.clock()

   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs  = self:rkStepperFields()[outIdx]
   
   self.calcPhi(tCurr, population, potCurr, potRhs)

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

-- Solve for dApardt in p>=2, or solve for a provisional dApardtProv in p=1.
function GkField:advanceStep2(tCurr, population, inIdx, outIdx)
   local tmStart = Time.clock()

   local potCurr = self:rkStepperFields()[inIdx]
   local potRhs  = self:rkStepperFields()[outIdx]

   if self.evolve then

      self.currentDensLocal:clear(0.0)
      self.modWeightAmpereLocal:copy(self.dApardtSlvr_lapModFac)
      for _, s in population.iterLocal() do
         if s:isEvolving() then 
            self.modWeightAmpereLocal:accumulate((s:getCharge()^2)/s:getMass(), s:getNumDensity())
            -- Taking momDensity at outIdx gives momentum moment of df/dt.
            self.currentDensLocal:accumulate(s:getCharge(), s:getMomDensity(outIdx))
         end
      end
      -- Reduce currents and Ampere weight across species communicator.
      self.currentDens:clear(0.0)
      population:AllreduceByCell(self.currentDensLocal, self.currentDens, 'sum')
      self.modWeightAmpere:clear(0.0)
      population:AllreduceByCell(self.modWeightAmpereLocal, self.modWeightAmpere, 'sum')

      self.dApardtSlvr:setModifierWeight(self.modWeightAmpere)
      -- dApar/dt solve.
      local dApardt = potCurr.dApardt
      self.dApardtSlvr:advance(tCurr, {self.currentDens}, {dApardt}) 

      -- Apply BCs.
      local tmStart = Time.clock()
      dApardt:sync(true)
      self.timers.bc = self.timers.bc + Time.clock() - tmStart

      -- Apar RHS is just dApar/dt.
      potRhs.apar:copy(dApardt)
   end

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

-- NOTE: global boundary conditions handled by solver. This just updates interproc ghosts.
-- Also NOTE: this method does not usually get called (because it is not called in applyBcIdx).
function GkField:applyBc(tCurr, potIn)
   local tmStart = Time.clock()

   potIn.phi:sync(true)
   if self.isElectromagnetic then 
      potIn.apar:sync(true) 
      potIn.dApardt:sync(true) 
   end

   self.timers.bc = self.timers.bc + Time.clock() - tmStart
end
   
function GkField:clearTimers() 
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
end

function GkField:printDevDiagnostics() 
   if self.phiSlvr then
      self.phiSlvr.printDevDiagnostics()
   end
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

   self.evolve = xsys.pickBool(tbl.evolve, false) -- by default these fields are not time-dependent.

   -- Create triggers to write fields.
   local nFrame = tbl.nFrame or appTbl.nFrame
   self.ioTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   self.ioFrame = 0 -- Frame number for IO.
   
   -- Get function to initialize background magnetic field.
   self.bmagFunc = tbl.bmag
   --assert(self.bmagFunc and type(self.bmagFunc)=="function", "GkGeometry: must specify background magnetic field function with 'bmag'")

   -- File containing geometry quantities that go into equations.
   self.fromFile = tbl.fromFile

   -- write ghost cells on boundaries of global domain (for BCs)
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   self.timers = {advance = 0.,   bc = 0.}
end

function GkGeometry:setGrid(grid)
   self.grid = grid
   self.ndim = self.grid:ndim()
end

function GkGeometry:alloc()
   -- Allocate fields.
   self.geo = {}

   local ghostNum     = {1,1}
   local syncPeriodic = false

   -- Background magnetic field.
   self.geo.bmag = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- bmagInv ~ 1/B.
   self.geo.bmagInv   = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   self.geo.bmagInvSq = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- cmag = J B / sqrt(g_zz).
   self.geo.cmag = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   -- Functions for Laplacian.
   -- g^xx = |grad x|**2.
   self.geo.gxx = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- g^xy = grad x . grad y.
   self.geo.gxy = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- g^yy = |grad y|**2.
   self.geo.gyy = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

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
   -- b_x = g_xz / ( sqrt(g_zz) )
   self.geo.b_x = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- b_y = g_yz / ( sqrt(g_zz) )
   self.geo.b_y = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   -- b_z = sqrt(g_zz)
   self.geo.b_z = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)

   self.geo.b_i = createField(self.grid,self.basis,ghostNum,3,syncPeriodic)
 
   -- Functions for laplacian, including Jacobian factor.
   self.geo.gxxJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   self.geo.gxyJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   self.geo.gyyJ = createField(self.grid,self.basis,ghostNum,1,syncPeriodic)
   
   if self.fromFile == nil then
      self.geo.allGeo = createField(self.grid,self.basis,ghostNum,15,syncPeriodic)
   end
end

function GkGeometry:createSolver(population)

   -- Get configuration-space Jacobian from grid (used by App/Projection/GkProjection, via GkSpecies).
   self.jacobGeoFunc = function (t, xn) return self.grid:calcJacobian(xn) end

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
         elseif self.grid._inDim == nil then -- Input file doesn't have mapc2p.
            g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
         end
         local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
         local b_x      = g_xz/math.sqrt(g_zz)
         local b_y      = g_yz/math.sqrt(g_zz)
         local b_z      = math.sqrt(g_zz)

         local det = jacobian^2
         local gxx = (g_yy*g_zz-g_yz^2)/det
         local gxy = (g_xz*g_yz-g_xy*g_zz)/det
         local gxz = (g_xy*g_yz-g_xz*g_yy)/det
         local gyy = (g_xx*g_zz-g_xz^2)/det
         local gyz = (g_xy*g_xz-g_xx*g_yz)/det
         local gzz = (g_xx*g_yy-g_xy^2)/det

         local bmag = self.bmagFunc(t, xn)
         local cmag = jacobian*bmag/math.sqrt(g_zz)

         return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, cmag, 
                b_x, b_y, b_z, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
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
         elseif self.grid._inDim == nil then -- Input file doesn't have mapc2p.
            g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
         end
         local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
         local b_x      = g_xz/math.sqrt(g_zz)
         local b_y      = g_yz/math.sqrt(g_zz)
         local b_z      = math.sqrt(g_zz)

         local det = jacobian^2
         local gxx = (g_yy*g_zz-g_yz^2)/det
         local gxy = (g_xz*g_yz-g_xy*g_zz)/det
         local gxz = (g_xy*g_yz-g_xz*g_yy)/det
         local gyy = (g_xx*g_zz-g_xz^2)/det
         local gyz = (g_xy*g_xz-g_xx*g_yz)/det
         local gzz = (g_xx*g_yy-g_xy^2)/det

         local bmag = self.bmagFunc(t, xn)
         local cmag = jacobian*bmag/math.sqrt(g_zz)

         return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, cmag, 
                b_x, b_y, b_z, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
       end
   else
      self.calcAllGeo = function(t, xn)
         local g = {}
         self.grid:calcMetric(xn, g)
         local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = g[1], g[2], g[3], g[4], g[5], g[6]
         local jacobian = math.sqrt(-g_xz^2*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz^2 - g_xy^2*g_zz + g_xx*g_yy*g_zz)
         local b_x      = g_xz/math.sqrt(g_zz)
         local b_y      = g_yz/math.sqrt(g_zz)
         local b_z      = math.sqrt(g_zz)

         local det = jacobian^2
         local gxx = (g_yy*g_zz-g_yz^2)/det
         local gxy = (g_xz*g_yz-g_xy*g_zz)/det
         local gxz = (g_xy*g_yz-g_xz*g_yy)/det
         local gyy = (g_xx*g_zz-g_xz^2)/det
         local gyz = (g_xy*g_xz-g_xx*g_yz)/det
         local gzz = (g_xx*g_yy-g_xy^2)/det

         local bmag = self.bmagFunc(t, xn)
         local cmag = jacobian*bmag/math.sqrt(g_zz)

         return jacobian, 1/jacobian, jacobian*bmag, 1/(jacobian*bmag), bmag, cmag, 
                b_x, b_y, b_z, gxx, gxy, gyy, gxx*jacobian, gxy*jacobian, gyy*jacobian
       end
   end

   if self.fromFile == nil then
      self.setAllGeo = Updater.EvalOnNodes {
         onGrid = self.grid,   evaluate = self.calcAllGeo,
         basis  = self.basis,  onGhosts = true,
      }
   end

   -- Updater to separate vector components packed into a single CartField.
   self.separateComponents = Updater.SeparateVectorComponents {
      onGrid = self.grid,  basis = self.basis,
   }
      
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.geo.bmag:elemType(),
      writeGhost = self.writeGhost,
      metaData   = { polyOrder = self.basis:polyOrder(),
	             basisType = self.basis:id(),
                     grid      = GKYL_OUT_PREFIX .. "_grid.bp", },
      writeRankInComm = {0, population:getComm_host(),}
   }   
end

function GkGeometry:createDiagnostics() end

function GkGeometry:initField(population)
   local log = Logger { logToFile = true }
   log("...Initializing the geometry...\n")
   if self.fromFile then
      -- Read the geometry quantities from a file.
      local tm, fr = self.fieldIo:read({jacobGeo=self.geo.jacobGeo, jacobGeoInv=self.geo.jacobGeoInv, jacobTot=self.geo.jacobTot,
         jacobTotInv=self.geo.jacobTotInv, bmag=self.geo.bmag,
         cmag=self.geo.cmag, b_x=self.geo.b_x, b_y=self.geo.b_y, b_z=self.geo.b_z, gxx=self.geo.gxx,
         gxy=self.geo.gxy, gyy=self.geo.gyy, gxxJ=self.geo.gxxJ, gxyJ=self.geo.gxyJ, gyyJ=self.geo.gyyJ},
         self.fromFile, true)
   else
      self.setAllGeo:advance(0.0, {}, {self.geo.allGeo})
      self.separateComponents:advance(0, {self.geo.allGeo},
         {self.geo.jacobGeo, self.geo.jacobGeoInv, self.geo.jacobTot, self.geo.jacobTotInv,
          self.geo.bmag, self.geo.cmag, self.geo.b_x, self.geo.b_y, self.geo.b_z,
          self.geo.gxx, self.geo.gxy, self.geo.gyy, self.geo.gxxJ, self.geo.gxyJ, self.geo.gyyJ})
   end
   local numB = self.basis:numBasis()
   self.geo.b_i:combineOffset(1, self.geo.b_x, 0*numB, 1, self.geo.b_y, 1*numB, 1, self.geo.b_z, 2*numB)

   -- Compute 1/B and 1/(B^2). LBO collisions require that this is
   -- done via weak operations instead of projecting 1/B or 1/B^2.
   local unitField = createField(self.grid,self.basis,{1,1})
   unitField:shiftc(1.*(math.sqrt(2.)^(self.ndim)), 0)
   local confWeakMultiply = Updater.CartFieldBinOp {
      weakBasis = self.basis,  operation = "Multiply",
      onGhosts  = false,
   }
   local confWeakDivide = Updater.CartFieldBinOp {
      weakBasis = self.basis,  operation = "Divide",
      onRange   = unitField:localRange(),  onGhosts = false,
   }
   confWeakDivide:advance(0., {self.geo.bmag, unitField}, {self.geo.bmagInv})
   confWeakMultiply:advance(0., {self.geo.bmagInv, self.geo.bmagInv}, {self.geo.bmagInvSq})

   log("...Finished initializing the geometry\n")

   -- Sync ghost cells. These calls do not enforce periodicity because
   -- these fields were created with syncPeriodicDirs = false.
   self.geo.bmag:sync(false)
   self.geo.bmagInv:sync(false)
   self.geo.bmagInvSq:sync(false)
   self.geo.cmag:sync(false)
   self.geo.gxx:sync(false)
   self.geo.gxy:sync(false)
   self.geo.gyy:sync(false)
   self.geo.gxxJ:sync(false)
   self.geo.gxyJ:sync(false)
   self.geo.gyyJ:sync(false)
   self.geo.jacobGeo:sync(false)
   self.geo.jacobGeoInv:sync(false)
   self.geo.jacobTotInv:sync(false)
   self.geo.b_x:sync(false)
   self.geo.b_y:sync(false)
   self.geo.b_z:sync(false)
   self.geo.b_i:sync(false)
   
   -- Apply BCs.
   self:applyBc(0, self.geo)
end

function GkGeometry:write(tm)
   -- Not evolving geometry, so only write geometry at beginning.
   if self.ioFrame == 0 then
      -- Write the geometry quantities to a file.
      for _, v in ipairs({{"%d",self.writeGhost},{"restart",true}}) do
         self.fieldIo:write({jacobGeo=self.geo.jacobGeo, jacobGeoInv=self.geo.jacobGeoInv, jacobTot=self.geo.jacobTot,
            jacobTotInv=self.geo.jacobTotInv, bmag=self.geo.bmag, bmagInv=self.geo.bmagInv,
            cmag=self.geo.cmag, b_x=self.geo.b_x, b_y=self.geo.b_y, b_z=self.geo.b_z, b_i=self.geo.b_i, gxx=self.geo.gxx,
            gxy=self.geo.gxy, gyy=self.geo.gyy, gxxJ=self.geo.gxxJ, gxyJ=self.geo.gxyJ, gyyJ=self.geo.gyyJ},
            string.format("allGeo_"..v[1]..".bp", self.ioFrame), tm, self.ioFrame, v[2])
      end

--      -- MF 2023/09/13: I don't remember why this was necessary if PlasmaOnCartGrid already calls the :write method of the grid.
--      -- Write a grid file.
--      local gridNodalCoords = self.grid:getNodalCoords()
--      if gridNodalCoords then self.fieldIo:write(gridNodalCoords, "grid.bp", tm, self.ioFrame) end
   end
   self.ioFrame = self.ioFrame+1
end

function GkGeometry:writeRestart(tm) end

function GkGeometry:rkStepperFields()
   return { self.geo, self.geo, self.geo, self.geo }
end

function GkGeometry:advance(tCurr) end

function GkGeometry:applyBcIdx(tCurr, idx)
   self:applyBc(tCurr, self:rkStepperFields()[1])
end

function GkGeometry:applyBc(tCurr, geoIn) end

function GkGeometry:clearTimers() 
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
end

return {GkField = GkField, GkGeometry = GkGeometry}
