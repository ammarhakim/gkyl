-- Gkyl ------------------------------------------------------------------------
--
-- App support code: GkProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto            = require "Lib.Proto"
local ProjectionBase   = require "App.Projection.ProjectionBase"
local xsys             = require "xsys"
local Updater          = require "Updater"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local DataStruct       = require "DataStruct"
local lume             = require "Lib.lume"

-- Shell class for gyrokinetic projections.
local GyrokineticProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GyrokineticProjection:init(tbl) self.tbl = tbl end

function GyrokineticProjection:fullInit(mySpecies)
   self.mass   = mySpecies:getMass()
   self.charge = mySpecies:getCharge()

   self.fromFile = self.tbl.fromFile

   self.exactScaleM0    = xsys.pickBool(self.tbl.exactScaleM0, true)
   self.exactScaleM012  = xsys.pickBool(self.tbl.exactScaleM012, false)
   if self.exactScaleM012 then self.exactScaleM0 = false end

   self.power = self.tbl.power
   self.scaleWithSourcePower = xsys.pickBool(self.tbl.scaleWithSourcePower, false)
end

function GyrokineticProjection:createSolver(mySpecies)
   self.phaseBasis = mySpecies.basis
   self.phaseGrid  = mySpecies.grid
   self.confBasis  = mySpecies.confBasis
   self.confGrid   = mySpecies.confGrid

   self.cdim = self.confGrid:ndim()
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.weakMultiplyConfPhase = Updater.CartFieldBinOp {
      weakBasis  = self.phaseBasis,  operation = "Multiply",
      fieldBasis = self.confBasis,   onGhosts  = true,
   }
end

--------------------------------------------------------------------------------
local FunctionProjection = Proto(GyrokineticProjection)

function FunctionProjection:createSolver(mySpecies)
   FunctionProjection.super.createSolver(self, mySpecies)

   self.initFunc = self.tbl.func
   if not self.initFunc then self.initFunc = self.tbl[1] end

   assert(self.initFunc, "FunctionProjection: Must specify the function")
   assert(type(self.initFunc) == "function", "The input must be a table containing function")

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType   = mySpecies.distf[1]:elemType(),
         method     = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData   = {polyOrder = self.phaseBasis:polyOrder(),
                       basisType = self.phaseBasis:id()}
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,   evaluate = self.initFunc,
         basis  = self.phaseBasis,  onGhosts = true
      }
   end
end

function FunctionProjection:allocConfField(metaData)
   local m = DataStruct.Field {
        onGrid        = self.confGrid,
        numComponents = self.confBasis:numBasis(),
        ghost         = {1, 1},
        metaData      = {polyOrder = self.confBasis:polyOrder(),
                         basisType = self.confBasis:id()},
   }
   m:clear(0.0)
   return m
end

function FunctionProjection:advance(time, inFlds, outFlds)
   local mySpecies, extField = inFlds[1], inFlds[2]
   local distf = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      -- Multiply by the phase jacobian. Unclear whether we should weak multiply or multiply
      -- within the projection. We choose the latter for now (since that's what's also done
      -- in MaxwellianProjection/MaxwellianOnBasis.
      local initFunc = self.initFunc
      if extField then
	 if extField.bmagFunc then
            initFunc = function(t, xn)
               local xconf = {}
               for d = 1, self.cdim do xconf[d] = xn[d] end
               local J = extField.bmagFunc(t,xconf) -- Phase space Jacobian.
               local f = self.initFunc(t,xn)
               return J*f
            end
         end
      end
      self.project:setFunc(initFunc)
      self.project:advance(time, {}, {distf})
   end
end

function FunctionProjection:scaleDensity(mySpecies, distf, currentM0, targetM0)
   local M0mod = mySpecies:allocMoment()

   local weakDivision = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,         operation = "Divide",
      onRange   = M0mod:localExtRange(),  onGhosts  = true,
   }

   -- Calculate M0mod = targetM0 / currentM0.
   weakDivision:advance(0.0, {currentM0, targetM0}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end

function FunctionProjection:createCouplingSolver(species, field, externalField)
   if not self.fromFile then
      local mySpecies = nil
      for _, s in lume.orderedIter(species) do
         if s.charge == self.charge then mySpecies = s end
      end

      if self.charge < 0.0 then
         -- Scale the electrons to have the same density as the ions.
         local numDens, numDensScaleTo = self:allocConfField(), self:allocConfField()
         local ionName, elcName = nil, nil
         for nm, _ in lume.orderedIter(species) do
            if 0. < s.charge then ionName = nm end
            if s.charge < 0. then elcName = nm end
         end
         species[elcName].numDensityCalc:advance(0.0, {species[elcName]:getDistF()}, {numDens})
         species[ionName].numDensityCalc:advance(0.0, {species[ionName]:getDistF()}, {numDensScaleTo})
         self:scaleDensity(species[elcName], species[elcName]:getDistF(), numDens, numDensScaleTo)
      end
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then self.weakMultiplyConfPhase:advance(0, {mySpecies:getDistF(), jacobGeo}, {mySpecies:getDistF()}) end
   end
end

--------------------------------------------------------------------------------
local MaxwellianProjection = Proto(GyrokineticProjection)

function MaxwellianProjection:createSolver(mySpecies)
   MaxwellianProjection.super.createSolver(self, mySpecies)

   local tbl = self.tbl
   self.density     = assert(tbl.density, "Maxwellian: must specify 'density'")
   self.driftSpeed  = tbl.driftSpeed or function(t, zn) return 0. end
   self.temperature = assert(tbl.temperature,
                             "Maxwellian: must specify 'temperature'")

   -- Check for constants instead of functions.
   if type(self.density) == "number" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.driftSpeed) == "number" then
      self.driftSpeed = function (t, zn) return tbl.driftSpeed end
   end
   if type(self.temperature) == "number" then
      self.temperature = function (t, zn) return tbl.temperature end
   end

   self.initFunc = function (t, zn)
      return mySpecies:Maxwellian(zn, self.density(t, zn),
                                      self.driftSpeed(t, zn),
                                      self.temperature(t, zn))
   end

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType  = mySpecies.distf[1]:elemType(),
         method    = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData  = {polyOrder = self.phaseBasis:polyOrder(),
                      basisType = self.phaseBasis:id(),
                      charge    = self.charge,
                      mass      = self.mass,},
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,   evaluate = self.initFunc,
         basis  = self.phaseBasis,  onGhosts = true
      }
   end
end

function MaxwellianProjection:allocConfField(vComp)
   local vComp = vComp or 1
   local m = DataStruct.Field {
        onGrid        = self.confGrid,
        numComponents = vComp*self.confBasis:numBasis(),
        ghost         = {1, 1},
        metaData      = {polyOrder = self.confBasis:polyOrder(),
                         basisType = self.confBasis:id()},
   }
   m:clear(0.0)
   return m
end

function MaxwellianProjection:scaleDensity(mySpecies, distf)
   local M0e, M0 = mySpecies:allocMoment(), mySpecies:allocMoment()
   local M0mod   = mySpecies:allocMoment()

   mySpecies.numDensityCalc:advance(0.0, {distf}, {M0})
   local project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,   evaluate = function (t, zn) return self.density(t, zn) end,
      basis  = self.confBasis,  onGhosts = true,
   }
   project:advance(0.0, {}, {M0e})

   local weakDivision = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Divide",
      onRange   = M0e:localExtRange(),  onGhosts  = true,
   }

   -- Calculate M0mod = M0e / M0.
   weakDivision:advance(0.0, {M0, M0e}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end

function MaxwellianProjection:scaleM012(mySpecies, distf)
   local sp                                        = mySpecies
   local M0, M2par, M2perp                         = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
   local M0_e, M2_e                                = sp:allocMoment(), sp:allocMoment()
   local M0_mod, M2par_mod, M2perp_mod             = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
   local M0par_mod, M0perp_mod                     = sp:allocMoment(), sp:allocMoment() 
   local M0par_mod2, M0perp_mod2                   = sp:allocMoment(), sp:allocMoment() 
   local M02par_mod, M02perp_mod                   = sp:allocMoment(), sp:allocMoment() 
   local distf0_mod, distf2par_mod, distf2perp_mod = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()

   -- Initialize maxwellian distribution distf0 = FM, along with 
   -- distf2par = vpar^2/2*FM and distf2perp = (mu*B/mass)*FM.
   local distf0, distf2par, distf2perp = sp:allocDistf(), sp:allocDistf(), sp:allocDistf()
   distf0:copy(distf)
   local phaseProject = Updater.ProjectOnBasis {
      onGrid   = self.phaseGrid,
      basis    = self.phaseBasis,
      evaluate = function(t,xn) return 0. end,   -- Set below.
      onGhosts = true
   }
   local distf2parFunc = function (t, zn)
      local xconf = {}
      for d = 1, self.cdim do xconf[d] = zn[d] end
      local vpar = zn[self.cdim+1]
      return vpar^2/2*sp.jacobPhaseFunc(t,xconf)*self.initFunc(t,zn)
   end
   phaseProject:setFunc(distf2parFunc)
   phaseProject:advance(0.0, {}, {distf2par})
   if self.vdim > 1 then 
      local distf2perpFunc = function (t, zn)
         local xconf = {}
         for d = 1, self.cdim do xconf[d] = zn[d] end
         local mu = zn[self.cdim+2]
         return mu*sp.bmagFunc(t,zn)/sp.mass*sp.jacobPhaseFunc(t,xconf)*self.initFunc(t,zn)
      end
      phaseProject:setFunc(distf2perpFunc)
      phaseProject:advance(0.0, {}, {distf2perp})
   end

   -- Calculate (inexact) moments of initial distribution function.
   sp.numDensityCalc:advance(0.0, {distf}, {M0})
   sp.M2parCalc:advance(0.0, {distf}, {M2par})
   if self.vdim > 1 then sp.M2perpCalc:advance(0.0, {distf}, {M2perp}) end

   -- Initialize exact moments.
   local confProject = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = function(t,xn) return 0. end,   -- Set below.
      onGhosts = true,
   }
   local M0func = function (t, zn)
      return self.density(t, zn)
   end
   confProject:setFunc(M0func)
   confProject:advance(0.0, {}, {M0_e})

   local M2func = function (t, zn)
      return self.density(t, zn)*self.temperature(t, zn, sp)/sp.mass
   end
   confProject:setFunc(M2func)
   confProject:advance(0.0, {}, {M2_e})

   -- Initialize weak multiplication/division operators.
   local weakDivision = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Divide",
      onRange   = M0_e:localExtRange(),  onGhosts  = true,
   }
   local weakMultiplicationConf = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
      onGhosts  = true,
   }
   local weakMultiplicationPhase = Updater.CartFieldBinOp {
      weakBasis  = self.phaseBasis,  operation  = "Multiply",
      fieldBasis = self.confBasis,   onGhosts   = true,
   }

   -- Calculate M0_mod = M0_e / M0.
   weakDivision:advance(0.0, {M0, M0_e}, {M0_mod})
   -- Calculate M2par_mod = M2_e / M2par.
   weakDivision:advance(0.0, {M2par, M2_e}, {M2par_mod})
   -- Calculate M2perp_mod = M2_e / M2perp.
   if self.vdim > 1 then weakDivision:advance(0.0, {M2perp, M2_e}, {M2perp_mod}) end
   -- Calculate M0par_mod = M0_e / M2par.
   weakDivision:advance(0.0, {M2par, M0_e}, {M0par_mod})
   -- Calculate M0par_mod2 = M0 / M2par.
   weakDivision:advance(0.0, {M2par, M0}, {M0par_mod2})
   if self.vdim > 1 then 
      -- Calculate M0perp_mod = M0_e / M2perp.
      weakDivision:advance(0.0, {M2perp, M0_e}, {M0perp_mod})
      -- Calculate M0perp_mod2 = M0 / M2perp.
      weakDivision:advance(0.0, {M2perp, M0}, {M0perp_mod2})
   end
   
   -- Calculate M02par_mod = M0par_mod2 * M2par_mod = (M0/M2par)*(M2_e/M2par).
   weakMultiplicationConf:advance(0.0, {M0par_mod2, M2par_mod}, {M02par_mod})
   -- Calculate M02perp_mod = M0perp_mod2 * M2perp_mod = (M0/M2perp)*(M2perp_e/M2perp).
   if self.vdim > 1 then weakMultiplicationConf:advance(0.0, {M0perp_mod2, M2perp_mod}, {M02perp_mod}) end

   -- Calculate distf modifiers from combinations of moment modifiers.
   if self.vdim==1 then 
      distf0_mod:combine(3/2, M0_mod, -1/2, M2par_mod)
   else 
      distf0_mod:combine(5/2, M0_mod, -1/2, M2par_mod, -1, M2perp_mod)
   end
   distf2par_mod:combine(1, M02par_mod, -1, M0par_mod)
   if self.vdim > 1 then distf2perp_mod:combine(1, M02perp_mod, -1, M0perp_mod) end

   -- Calculate distf0 = distf0_mod * distf0.
   weakMultiplicationPhase:advance(0.0, {distf0_mod, distf0}, {distf0})
   -- Calculate distf2par = distf2par_mod * distf2par.
   weakMultiplicationPhase:advance(0.0, {distf2par_mod, distf2par}, {distf2par})
   -- Calculate distf2perp = distf2perp_mod * distf2perp.
   if self.vdim > 1 then weakMultiplicationPhase:advance(0.0, {distf2perp_mod, distf2perp}, {distf2perp}) end

   -- Combine and finish.
   distf:combine(1, distf0, 1, distf2par)
   if self.vdim > 1 then distf:accumulate(1, distf2perp) end
end

function MaxwellianProjection:advance(time, inFlds, outFlds)
   local mySpecies, extField = inFlds[1], inFlds[2]
   local distf = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      local bmag = extField.geo.bmag
      -- Project the moments onto configuration-space basis.
      local confScalarProject = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 0. end,   -- Set below.
         basis  = self.confBasis,  onGhosts = true
      }
      local confVec2Project = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 0., 0. end,   -- Set below.
         basis  = self.confBasis,  onGhosts = true
      }
      local numDens  = self:allocConfField()
      local primMoms = self:allocConfField(2)
      if self.exactScaleM0 then
         -- Use a unit density because we are going to rescale the density anyways,
         -- and it is easier to weak-divide by something close to unity.
         confScalarProject:setFunc(function(t, xn) return 1. end)
      else
         confScalarProject:setFunc(self.density)
      end
      confScalarProject:advance(time, {}, {numDens})
      confVec2Project:setFunc(function(t,xn) return self.driftSpeed(t,xn), self.temperature(t,xn)/self.mass end)
      confVec2Project:advance(time, {}, {primMoms})
      -- Project the Maxwellian. It includes a factor of jacobPhase=B*_||.
      local projMaxwell = Updater.MaxwellianOnBasis {
         onGrid      = self.phaseGrid,   confBasis  = self.confBasis,
         phaseBasis  = self.phaseBasis,  mass       = self.mass,
         usePrimMoms = true,             onGhosts   = true,
      }
      -- Use bmag as the total jacobian here because we weak multiply by jacobGeo later.
      projMaxwell:advance(time,{numDens,primMoms,bmag,bmag},{distf})

      if self.exactScaleM0 then
         self:scaleDensity(mySpecies, distf)
      elseif self.exactScaleM012 then
         self:scaleM012(mySpecies, distf)
      end

      local jacobGeo = extField.geo.jacobGeo
      if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
   end
end


----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
