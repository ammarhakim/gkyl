-- Gkyl ------------------------------------------------------------------------
--
-- App support code: GkProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local Updater    = require "Updater"
local xsys       = require "xsys"
local DataStruct = require "DataStruct"
local FunctionProjectionParent   = require ("App.Projection.KineticProjection").FunctionProjection
local MaxwellianProjectionParent = require ("App.Projection.KineticProjection").MaxwellianProjection

--------------------------------------------------------------------------------
-- Gk-specific GkProjection.FunctionProjection includes Jacobian factors in initFunc.
local FunctionProjection = Proto(FunctionProjectionParent)

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
   local extField = inFlds[1]
   local distf    = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
--      if self.species.jacobPhaseFunc and self.vdim > 1 then
--         local initFuncWithoutJacobian = self.initFunc
--         self.initFunc = function (t, xn)
--            local xconf = {}
--            for d = 1, self.cdim do xconf[d] = xn[d] end
--            local J = self.species.jacobPhaseFunc(t,xconf)
--            local f = initFuncWithoutJacobian(t,xn)
--            return J*f
--         end
--      end

      -- Note: don't use self.project as this does not have jacobian factors in initFunc.
      local project = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,   evaluate = self.initFunc,
         basis  = self.phaseBasis,  onGhosts = true
      }
      project:advance(time, {}, {distf})
   end

   local jacobTot, jacobPhase = extField.geo.jacobTot, extField.geo.bmag
   if jacobTot then self.weakMultiplyConfPhase:advance(0, {distf, jacobTot}, {distf})
   elseif jacobPhase then self.weakMultiplyConfPhase:advance(0, {distf, jacobPhase}, {distf}) end
end

function FunctionProjection:createCouplingSolver(species,field, externalField)
   if not self.fromFile then
      if self.species.charge < 0.0 then
         -- Scale the electrons to have the same density as the ions.
         local numDens = self:allocConfField()
         local numDensScaleTo = self:allocConfField()
         local ionName = nil
         for nm, s in lume.orderedIter(species) do
            if 0.0 < s.charge then ionName = nm end
         end
         self.species.numDensityCalc:advance(0.0, {self.species:getDistF()}, {numDens})
         species[ionName].numDensityCalc:advance(0.0, {species[ionName]:getDistF()}, {numDensScaleTo})
         self:scaleDensity(self.species:getDistF(), numDens, numDensScaleTo)
      end
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then self.weakMultiplyConfPhase:advance(0, {self.species:getDistF(), jacobGeo}, {self.species:getDistF()}) end
   end
end

--------------------------------------------------------------------------------
-- Gk-specific GkProjection.MaxwellianProjection extends MaxwellianProjection base class, including 
-- adding jacobian factors in initFunc.
local MaxwellianProjection = Proto(MaxwellianProjectionParent)

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

function MaxwellianProjection:lagrangeFix(distf)
   local M0, dM0 = self.species:allocMoment(), self.species:allocMoment()
   local M1, dM1 = self.species:allocMoment(), self.species:allocMoment()
   local M2, dM2 = self.species:allocMoment(), self.species:allocMoment()

   local project = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = function(t,xn) return 0. end,   -- Set below.
      onGhosts = true,
   }

   self.species.numDensityCalc:advance(0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   project:setFunc(func)
   project:advance(0.0, {}, {dM0})
   dM0:accumulate(-1.0, M0)

   self.species.momDensityCalc:advance(0.0, {distf}, {M1})
   func = function (t, zn)
      return self.density(t, zn, self.species) *
	 self.driftSpeed(t, zn, self.species)
   end
   project:setFunc(func)
   project:advance(0.0, {}, {dM1})
   dM1:accumulate(-1.0, M1)

   self.species.ptclEnergyCalc:advance(0.0, {distf}, {M2})
   func = function (t, zn)
      local drifts = self.driftSpeed(t, zn, self.species)
      if self.vdim == 1 then
	 return self.density(t, zn, self.species) *
	    (self.driftSpeed(t, zn, self.species)*self.driftSpeed(t, zn, self.species) + self.temperature(t, zn, self.species)/self.species.mass )
      else
	 return self.density(t, zn, self.species) *
	    (self.driftSpeed(t, zn, self.species)*self.driftSpeed(t, zn, self.species) + 3*self.temperature(t, zn, self.species)/self.species.mass )
      end
   end
   project:setFunc(func)
   project:advance(0.0, {}, {dM2})
   dM2:accumulate(-1.0, M2)

   local lagFix = Updater.LagrangeFix {
      onGrid     = self.phaseGrid,   confBasis  = self.confBasis,
      phaseBasis = self.phaseBasis,  mode       = 'gk',
      confGrid   = self.confGrid,    mass       = self.species.mass,
   }
   lagFix:advance(0.0, {dM0, dM1, dM2, self.species.bmag}, {distf})
end

function MaxwellianProjection:scaleM012(distf)
   local sp                                        = self.species
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
      return self.density(t, zn, sp)
   end
   confProject:setFunc(M0func)
   confProject:advance(0.0, {}, {M0_e})

   local M2func = function (t, zn)
      return self.density(t, zn, sp)*self.temperature(t, zn, sp)/sp.mass
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
   local extField = inFlds[1]
   local distf    = outFlds[1]
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
   end

   if self.exactScaleM0 then
      self:scaleDensity(distf)
   elseif self.exactScaleM012 then
      self:scaleM012(distf)
   end
   if self.exactLagFixM012 then self:lagrangeFix(distf) end

   local jacobGeo = extField.geo.jacobGeo
   if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
end


----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
