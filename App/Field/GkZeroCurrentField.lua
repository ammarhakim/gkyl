-- Gkyl ------------------------------------------------------------------------
--
-- Field solver from adding the q*vpar moments of the kinetic equation.
-- The result is an equation for Epar. Valid for quasinuetral plasmas only.
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
-- In order to compute dBdz.
local math = require("sci.math").generic
local diff = require("sci.diff")

local GkZeroCurrentField = Proto(FieldBase.FieldBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkZeroCurrentField:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkZeroCurrentField:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.ioMethod = "MPI"

   if appTbl.periodicDirs then self.periodicDirs = appTbl.periodicDirs else self.periodicDirs = {} end

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

   -- For storing integrated energies.
   self.EparSq = DataStruct.DynVector { numComponents = 1 }

   self.adiabatic = false  -- Electrons must be adiabatic. We check this later.

   -- Flag to indicate if Epar has been calculated.
   self.calcedEpar = false

   self._first = true

   self.bcTime = 0.0 -- Timer for BCs.
   self.timers = {advTime={0.,0.,0.}}
end

function GkZeroCurrentField:hasEB() return true, false end
function GkZeroCurrentField:setGrid(grid) self.grid = grid; self.ndim = self.grid:ndim() end

local function createField(grid, basis, ghostCells, vComp, periodicSync)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid           = grid,
      numComponents    = basis:numBasis()*vComp,
      ghost            = ghostCells,
      metaData         = {polyOrder = basis:polyOrder(),
                          basisType = basis:id()},
      syncPeriodicDirs = periodicSync,
   }
   fld:clear(0.0)
   return fld
end

function GkZeroCurrentField:alloc(nRkDup)
   self.Epar = createField(self.grid,self.basis,{1,1})
   self.EparAux = createField(self.grid,self.basis,{1,1})

   -- Allocate fields needed in RK update.
   -- nField is related to number of RK stages.
   self.fields = {}
   for i = 1, nRkDup do
      self.fields[i] = {}
      self.fields[i].Epar = self.Epar
      self.fields[i].EparAux = self.EparAux
   end

   -- Create fields for total charge density and div{P}.
   self.chargeDens = createField(self.grid,self.basis,{1,1})
   self.divP       = createField(self.grid,self.basis,{1,1})
   -- Fields to store M2par and M2perp moments of a distribution function.
   self.M2par  = createField(self.grid,self.basis,{1,1})
   self.M2perp = createField(self.grid,self.basis,{1,1})

end

-- Solve for initial fields self-consistently
-- from initial distribution function.
function GkZeroCurrentField:initField(species)
   -- Solve for initial phi.
   self:advance(0.0, species, 1, 1)
   self:phiSolve(0.0, species, 1, 1)

   -- Apply BCs and update ghosts.
   self:applyBc(0, self.fields[1])
end

function GkZeroCurrentField:rkStepperFields() return self.fields end

-- For RK timestepping for non-elliptic fields (e.g. only apar).
function GkZeroCurrentField:copyRk(outIdx, aIdx) end
function GkZeroCurrentField:combineRk(outIdx, a, aIdx, ...) end

function GkZeroCurrentField:createSolver(species, externalField)
   self.isElliptic = true   -- Needed so field is calculated self-consistently at end of full RK timestep.

   -- Get adiabatic species info.
   for nm, s in lume.orderedIter(species) do
      if Species.AdiabaticSpecies.is(s) then
         self.adiabatic, self.adiabSpec = true, s
      end
      -- Save ion and electron species names to be used later
      -- in re-scaling ions to enforce quasineutrality.
      if s.charge > 0. then self.ionName = nm elseif s.charge < 0. then self.elcName = nm end
   end

   self.recomputeElcPrimMom = false
   if species[self.elcName].needSelfPrimMom then
      self.recomputeElcPrimMom = true
      species[self.elcName].needSelfPrimMom = false -- So that electron self-prim mom are only computed once.
   end

   self.EparDivPSlvr = Updater.ZeroCurrentGkEparDivP {
      onGrid = self.grid,
      basis  = self.basis,
   }
   self.weakDiv = Updater.CartFieldBinOp {
      onGrid    = self.grid,   operation = "Divide",
      weakBasis = self.basis,  onGhosts  = false,
   }

   -- Set up constant dummy field.
   self.unitFld = createField(self.grid,self.basis,{1,1})
   local initUnit = Updater.ProjectOnBasis {
      onGrid = self.grid,   evaluate = function (t,xn) return 1.0 end,
      basis  = self.basis,  onGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitFld})

   -- We will need the varios geometrical quantities. Save pointers to them.
   local cmag, jacobTotInv, bmagInv
   if externalField.geo then
      if externalField.geo.name=="SimpleHelical" then
         cmag             = self.unitFld
         jacobTotInv      = self.unitFld
         self.jacobGeoInv = self.unitFld
      elseif externalField.geo.name=="GenGeo" then
         cmag             = externalField.geo.cmag
         jacobTotInv      = externalField.geo.jacobTotInv
         self.jacobGeoInv = externalField.geo.jacobGeoInv
      end
      bmagInv = externalField.geo.bmagInv
   else
      cmag             = self.unitFld
      jacobTotInv      = self.unitFld
      self.jacobGeoInv = self.unitFld
      bmagInv          = self.unitFld
   end

   -- Compute dBdz. Eventually we'll remove this or put it in geo.
   local dBdzFunc = function (t, xn)
      local function bmagUnpack(...)
         local xn1 = {...}
         return externalField.bmagFunc(0, xn1)
      end
      local deriv   = diff.derivativef(bmagUnpack, #xn)
      local xntable = {}
      for i = 1, #xn do xntable[i] = xn[i] end
      local f, dx = deriv(unpack(xntable))
      return dx
   end
   local dBdz = createField(self.grid,self.basis,{1,1})
   dBdz:clear(0.)
   local evOnNodes = Updater.EvalOnNodes {
      onGrid = self.grid,   evaluate = dBdzFunc,
      basis  = self.basis,  onGhosts = true,
   }
   evOnNodes:advance(0., {}, {dBdz})
   -- Compute d(ln B)/dz=(1/B)*dB/dz.
   self.weakMult = Updater.CartFieldBinOp {
      onGrid    = self.grid,   operation = "Multiply",
      weakBasis = self.basis,  onGhosts  = true,
   }
   self.dlnbmagdz = createField(self.grid,self.basis,{1,1})
   self.weakMult:advance(0., {bmagInv,dBdz}, {self.dlnbmagdz})

   -- Factor appearing in front of delpar.
   self.delparFac = createField(self.grid,self.basis,{1,1})
   self.weakMult:advance(0., {cmag,jacobTotInv}, {self.delparFac})

   -- The simulation appears to be unstable if we don't smooth Epar.
   -- MF: I expected discontinuous Epar to be ok since it only appears
   --     in the vpar surface term, but maybe I missed something.
   self.EparZSmoother = Updater.FemParPoisson {
      onGrid = self.grid,   bcLower = {{T="N",V=0.0}},
      basis  = self.basis,  bcUpper = {{T="N",V=0.0}},
      smooth = true,
   }

   self.denRat = createField(self.grid,self.basis,{1,1})

   self.nonPeriodicBCs = {}   -- List of non-periodic BCs to apply.

   -- Function to construct a BC updater.
   local function makeOpenBcUpdater(dir, edge)
      local bcOpen = function(dir, tm, idxIn, fIn, fOut)
         -- Requires skinLoop = "pointwise".
         self.basis:flipSign(dir, fIn, fOut)
      end

      return Updater.Bc {
         onGrid = self.grid,  edge = edge,
         dir    = dir,        boundaryConditions = {bcOpen},
         skinLoop = "pointwise",
      }
   end

   -- For non-periodic dirs, use open BCs. It's the most sensible choice given that the
   -- coordinate mapping could diverge outside of the interior domain.
   for dir = 1, self.ndim do
      if not lume.any(self.periodicDirs, function(t) return t==dir end) then
         self.nonPeriodicBCs["lower"] = makeOpenBcUpdater(dir, "lower")
         self.nonPeriodicBCs["upper"] = makeOpenBcUpdater(dir, "upper")
      end
   end

end

function GkZeroCurrentField:createDiagnostics()
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType   = self.fields[1].Epar:elemType(),
      method     = self.ioMethod,
      writeGhost = self.writeGhost,
      metaData   = {polyOrder = self.basis:polyOrder(),
                    basisType = self.basis:id(),},
   }
   -- Updaters for computing integrated quantities.
   self.intSqCalc = Updater.CartFieldIntegratedQuantCalc {
      onGrid   = self.grid,
      basis    = self.basis,
      quantity = "V2"
   }
end

function GkZeroCurrentField:write(tm, force)
   if self.calcIntFieldEnergyTrigger(tm) then
      -- Compute quantities integrated over domain.
      self.intSqCalc:advance(tm, { self.fields[1].Epar }, { self.EparSq })
   end

   if self.ioTrigger(tm) or force then
      self.fieldIo:write(self.fields[1].Epar, string.format("Epar_%d.bp", self.ioFrame), tm, self.ioFrame)
      self.EparSq:write(string.format("EparSq.bp"), tm, self.ioFrame)
 
      self.ioFrame = self.ioFrame+1
   end
end

function GkZeroCurrentField:writeRestart(tm)
   -- (the final "false" prevents writing of ghost cells).
   self.fieldIo:write(self.fields[1].Epar, "Epar_restart.bp", tm, self.ioFrame, false)

   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self.EparSq:write("EparSq_restart.bp", tm, self.dynVecRestartFrame, false, false)
   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function GkZeroCurrentField:readRestart()
   -- This read of restart file of Epar is only to get frame numbering correct.
   -- The forward Euler recomputes the potential before updating the hyperbolic part.
   local tm, fr = self.fieldIo:read(self.fields[1].Epar, "Epar_restart.bp")
   self.EparSq:read("EparSq_restart.bp", tm)

   self:applyBc(0, self.fields[1])

   self.ioFrame = fr
   -- Iterate triggers.
   self.ioTrigger(tm)
end

-- Solve for electrostatic potential Epar.
function GkZeroCurrentField:advance(tCurr, species, inIdx, outIdx)
   local tmStart = Time.clock()

   local fieldsCurr = self:rkStepperFields()[inIdx]
   local fieldsRhs  = self:rkStepperFields()[outIdx]

   -- ............................................................... --
   -- Rescale the electron density to match the ion density:
   self.weakDiv:advance(tCurr, {species[self.elcName]:getNumDensity(),species[self.ionName]:getNumDensity()}, {self.denRat}) 
   species[self.elcName].confPhaseWeakMultiply:advance(tCurr, {species[self.elcName]:rkStepperFields()[inIdx],self.denRat},
                                                       {species[self.elcName]:rkStepperFields()[inIdx]})
   -- Recompute coupling and primitive moments:
   if self.recomputeElcPrimMom then
      species[self.elcName].threeMomentsLBOCalc:advance(tCurr, {species[self.elcName]:rkStepperFields()[inIdx]}, { species[self.elcName].numDensity, species[self.elcName].momDensity, species[self.elcName].ptclEnergy,
                                                 species[self.elcName].m1Correction, species[self.elcName].m2Correction,
                                                 species[self.elcName].m0Star, species[self.elcName].m1Star, species[self.elcName].m2Star })
      -- Also compute species[self.elcName].primitive moments uPar and vtSq.
      species[self.elcName].primMomSelf:advance(tCurr, {species[self.elcName].numDensity, species[self.elcName].momDensity, species[self.elcName].ptclEnergy,
                                         species[self.elcName].m1Correction, species[self.elcName].m2Correction,
                                         species[self.elcName].m0Star, species[self.elcName].m1Star, species[self.elcName].m2Star}, {species[self.elcName].uParSelf, species[self.elcName].vtSqSelf})
      -- Indicate that moments, boundary corrections, star moments
      -- and self-primitive moments have been computed.
      for iF=1,4 do species[self.elcName].momentFlags[iF] = true end
   else
      species[self.elcName].numDensityCalc:advance(tCurr, {species[self.elcName]:rkStepperFields()[inIdx]}, { species[self.elcName].numDensity })
      -- Indicate that first moment has been computed.
      species[self.elcName].momentFlags[1] = true
   end
   -- ............. Finished rescaling to enforce ni=ne ............. --
   
   self.chargeDens:clear(0.0)
   self.divP:clear(0.0)
   for _, s in lume.orderedIter(species) do
      -- Add contribution to the total (q^2/m) weighted charge density.
      self.chargeDens:accumulate((s:getCharge()^2)/s:getMass(), s:getNumDensity())
      -- Compute the M2par and M2perp moments of this species.
      s.M2parCalc:advance(tCurr, {s:rkStepperFields()[inIdx]}, {self.M2par})
      s.M2perpCalc:advance(tCurr, {s:rkStepperFields()[inIdx]}, {self.M2perp})
      -- Add contribution to the divergence of the total pressure.
      self.EparDivPSlvr:advance(tCurr,{s:getCharge(),self.delparFac,self.dlnbmagdz,
                                       self.M2par,self.M2perp},{self.divP})
   end
   -- Divide hat{b} . div{ P } by delparFac*sum_s(q_s^2*J*n_s/m_s) 
   self.weakMult:advance(tCurr, {self.chargeDens,self.delparFac}, {self.chargeDens})
   self.weakDiv:advance(tCurr, {self.chargeDens,self.divP}, {fieldsCurr.EparAux})

   self.EparZSmoother:advance(tCurr, {fieldsCurr.EparAux}, {fieldsCurr.Epar})

   for _, bc in pairs(self.nonPeriodicBCs) do bc:advance(tCurr, {}, {fieldsCurr.Epar}) end
   fieldsCurr.Epar:sync(true)

   self.timers.advTime[1] = self.timers.advTime[1] + Time.clock() - tmStart
end

function GkZeroCurrentField:phiSolve(tCurr, species, inIdx, outIdx)
   -- Finish computing phi and apply BCs to phi.
   -- This could also be done in :advance.
end

-- NOTE: global boundary conditions handled by solver. This just updates interproc ghosts.
-- Also NOTE: this method does not usually get called (because it is not called in applyBcIdx).
function GkZeroCurrentField:applyBc(tCurr, fieldsIn)
   local tmStart = Time.clock()
   fieldsIn.Epar:sync(true)
   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end

function GkZeroCurrentField:totalSolverTime()
   local time = 0.
   if self.EparDivPSlvr then
      time = time + self.timers.advTime[1]
   end
   return time
end
function GkZeroCurrentField:totalBcTime() return self.bcTime end
function GkZeroCurrentField:energyCalcTime()
   local t = self.intSqCalc.totalTime
   return t
end

function GkZeroCurrentField:printDevDiagnostics() self.EparDivPSlvr:printDevDiagnostics() end

return GkZeroCurrentField
