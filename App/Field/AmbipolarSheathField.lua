-- Gkyl ------------------------------------------------------------------------
--
-- Field solver for sheath simulations with an ambipolar (logical) sheath
-- BC and adiabatic electrons.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local FieldBase        = require "App.Field.FieldBase"
local DataStruct       = require "DataStruct"
local LinearTrigger    = require "LinearTrigger"
local Proto            = require "Lib.Proto"
local Updater          = require "Updater"
local xsys             = require "xsys"
local FieldBase        = require "App.Field.FieldBase"
local Species          = require "App.Species"
local Time             = require "Lib.Time"
local lume             = require "Lib.lume"

local AmbipolarSheathField = Proto(FieldBase.FieldBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function AmbipolarSheathField:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function AmbipolarSheathField:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   -- Create triggers to write diagnostics.
   local nFrame = tbl.nFrame or appTbl.nFrame
   self.ioTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   -- Create trigger for how frequently to compute integrated diagnostics.
   -- Do not compute the integrated diagnostics less frequently than we output data.
   if appTbl.calcIntQuantEvery then
      self.calcIntFieldEnergyTrigger = LinearTrigger(0, appTbl.tEnd,  math.max(nFrame,math.floor(1/appTbl.calcIntQuantEvery)))
   else
      self.calcIntFieldEnergyTrigger = function(t) return true end
   end

   self.ioFrame = 0 -- Frame number for IO.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).

   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   self.adiabatic = false  -- Electrons must be adiabatic. We check this later.
   self.discontinuousPhi = xsys.pickBool(tbl.discontinuousPhi, false)

   self.bcTime = 0.0 -- Timer for BCs.

   self.timers = {advTime={0.,0.,0.}}
end

function AmbipolarSheathField:hasEB() return true, false end

function AmbipolarSheathField:setGrid(grid)
   self.grid = grid
   self.ndim = self.grid:ndim()
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

function AmbipolarSheathField:alloc(nRkDup)
   self.phi    = createField(self.grid,self.basis,{1,1})
   self.phiAux = createField(self.grid,self.basis,{1,1})

   self.zeroField = createField(self.grid,self.basis,{1,1})

   -- Allocate fields needed in RK update.
   -- nField is related to number of RK stages.
   self.potentials = {}
   for i = 1, nRkDup do
      self.potentials[i] = {}
      self.potentials[i].phi    = self.phi
      self.potentials[i].phiAux = self.phiAux
      self.potentials[i].apar    = self.zeroField
      self.potentials[i].dApardt = self.zeroField
   end

   -- For storing integrated energies.
   self.intPhiSq = DataStruct.DynVector { numComponents = 1, }
end

-- Solve for initial fields self-consistently
-- from initial distribution function.
function AmbipolarSheathField:initField(population)
   -- Solve for initial phi.
   self:advance(0.0, population, 1, 1)

   -- Apply BCs and update ghosts.
   self:applyBc(0, self.potentials[1])
end

function AmbipolarSheathField:rkStepperFields() return self.potentials end

-- For RK timestepping for non-elliptic fields (e.g. only apar).
function AmbipolarSheathField:copyRk(outIdx, aIdx) end
function AmbipolarSheathField:combineRk(outIdx, a, aIdx, ...) end

function AmbipolarSheathField:createSolver(population, externalField)
   self.isElliptic = true   -- Needed so field is calculated self-consistently at end of full RK timestep.
   
   local species = population:getSpecies()

   -- Get adiabatic species info.
   for nm, s in lume.orderedIter(species) do
      if Species.AdiabaticSpecies.is(s) then
         self.adiabatic, self.adiabSpec = true, s
      end
      if s.charge > 0. then self.ionName = nm else self.elcName = nm end
   end
   assert(self.adiabatic and (self.adiabSpec.charge<0.), "AmbipolarSheathField: currently only available for the adiabatic electron case.")

   local ndim = self.grid:ndim()
   local ionSpec = species[self.ionName]
   self.ionBC = {}
   local hasNonPeriodic = false 
   for _, bc in lume.orderedIter(ionSpec.nonPeriodicBCs) do
      if ionSpec.hasSheathBCs or bc.bcKind == "absorb" then
         if (ndim==3 and bc:getDir()==3) or ndim==1 then
            self.ionBC[bc:getEdge()] = bc
            self.ionBC[bc:getEdge()]:setSaveFlux(true)
            hasNonPeriodic = true
         end
      end
   end
   assert(hasNonPeriodic, "App.Field.AmbipolarSheathField: has to have sheath BCs.")
   lume.setOrder(self.ionBC)

   self.bcIonM0flux = {}
   for _, bc in lume.orderedIter(self.ionBC) do
      self.bcIonM0flux[bc:getEdge()] = bc:allocMoment()
   end

   self.phiSlvr = Updater.ASheathPotential {
      onGrid         = self.grid,  basis = self.basis,
      electronMass   = species[self.elcName]:getMass(),
      electronCharge = species[self.elcName]:getCharge(),
      electronTemp   = species[self.elcName]:temp(),
      boundaryGrids  = {lower=self.ionBC["lower"].confBoundaryGrid,
                        upper=self.ionBC["upper"].confBoundaryGrid},
   }
   self.phiZSmoother = Updater.FemParproj {
      onGrid              = self.grid,  basis   = self.basis, 
      periodicParallelDir = false,      onField = self:rkStepperFields()[1].phi,
   }

   -- We will need the reciprocal of the conf-space Jacobian (1/J) later.
   self.jacobGeoInv = externalField.geo.jacobGeoInv

   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.potentials[1].phi:elemType(),
      writeGhost = self.writeGhost,
      metaData   = {polyOrder = self.basis:polyOrder(),
                    basisType = self.basis:id(),},
      writeRankInComm = {0, population:getComm_host(),}
   }
end

function AmbipolarSheathField:createDiagnostics()
   -- Updaters for computing integrated quantities.
   self.intSqCalc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,   quantity = "V2",
      basis  = self.basis,
   }
end

function AmbipolarSheathField:write(tm, force)
   if self.calcIntFieldEnergyTrigger(tm) then
      -- Compute quantities integrated over domain.
      self.intSqCalc:advance(tm, { self.potentials[1].phi }, { self.intPhiSq })
   end

   if self.ioTrigger(tm) or force then
      self.fieldIo:write(self.potentials[1].phi, string.format("phi_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.intPhiSq:write(string.format("intPhiSq.bp"), tm, self.ioFrame)
 
      self.ioFrame = self.ioFrame+1
   end
end

function AmbipolarSheathField:writeRestart(tm)
   -- (the final "false" prevents writing of ghost cells).
   self.fieldIo:write(self.potentials[1].phi, "phi_restart.bp", tm, self.ioFrame, false)

   -- Write boundary fluxes. They are used upon restart.
   for _, bc in lume.orderedIter(self.ionBC) do
      bc:getBoundaryFluxFields()[1]:write(string.format("field_%s_restart.bp", bc.diagnostics.name))
   end

   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function AmbipolarSheathField:readRestart()
   -- This read of restart file of phi is only to get frame numbering correct.
   -- The forward Euler recomputes the potential before updating the hyperbolic part.
   local tm, fr = self.fieldIo:read(self.potentials[1].phi, "phi_restart.bp")

   self:applyBc(0, self.potentials[1])

   -- Read boundary particle fluxes from previous simulation.
   for _, bc in lume.orderedIter(self.ionBC) do
      local _, _ = self.fieldIo:read(bc:getBoundaryFluxFields()[1], string.format("field_%s_restart.bp", bc.diagnostics.name))
   end

   self.ioFrame = fr
   -- Iterate triggers.
   self.ioTrigger(tm)
end

-- Solve for electrostatic potential phi.
function AmbipolarSheathField:advance(tCurr, population, inIdx, outIdx)
   local tmStart = Time.clock()

   local potCurr = self:rkStepperFields()[inIdx]
   self.phiSlvr:advance(tCurr, {self.bcIonM0flux, population:getSpecies()[self.ionName]:getNumDensity(), self.jacobGeoInv},
                        {potCurr.phiAux})
   -- Smooth phi in z to ensure continuity in all directions.
   if self.ndim ~= 2 and not self.discontinuousPhi then
      self.phiZSmoother:advance(tCurr, {potCurr.phiAux}, {potCurr.phi})
   else
      potCurr.phi = potCurr.phiAux
   end

   -- Apply BCs. Make sure phi is continuous across skin-ghost boundary.
   local tmStart = Time.clock()
   potCurr.phi:sync(true)
   self.bcTime = self.bcTime + (Time.clock()-tmStart)

   self.timers.advTime[1] = self.timers.advTime[1] + Time.clock() - tmStart
end

function AmbipolarSheathField:useBoundaryFlux(tCurr, outIdx)
   -- Immediately after boundary fluxes are stored, use them to compute the particle
   -- flux to the boundary. Otherwise they get altered by :combineRk and multiplied by dt.
   for _, bc in lume.orderedIter(self.ionBC) do
      bc.numDensityCalc:advance(tCurr, {bc:getBoundaryFluxFields()[outIdx]}, {self.bcIonM0flux[bc:getEdge()]})
   end
end

-- NOTE: global boundary conditions handled by solver. This just updates interproc ghosts.
-- Also NOTE: this method does not usually get called (because it is not called in applyBcIdx).
function AmbipolarSheathField:applyBc(tCurr, potIn)
   local tmStart = Time.clock()
   potIn.phi:sync(true)
   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end

function AmbipolarSheathField:totalSolverTime()
   local time = 0.
   if self.phiSlvr then
      time = time + self.timers.advTime[1]+self.phiSlvr.totalTime
   end
   return time
end
function AmbipolarSheathField:totalBcTime() return self.bcTime end
function AmbipolarSheathField:energyCalcTime()
   local t = self.intSqCalc.totalTime
   return t
end

function AmbipolarSheathField:printDevDiagnostics() self.phiSlvr:printDevDiagnostics() end

return AmbipolarSheathField
