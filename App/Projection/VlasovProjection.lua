-- Gkyl ------------------------------------------------------------------------
--
-- App support code: VlasovProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

--local Time = require "Lib.Time".
local FunctionProjectionParent   = require ("App.Projection.KineticProjection").FunctionProjection
local MaxwellianProjectionParent = require ("App.Projection.KineticProjection").MaxwellianProjection
local Proto                      = require "Lib.Proto"
local Updater                    = require "Updater"
local xsys                       = require "xsys"

----------------------------------------------------------------------
-- Vlasov-specific VlasovProjection.FunctionProjection includes Jacobian factors ONLY for neutrals in general geometry.
local FunctionProjection = Proto(FunctionProjectionParent)

function FunctionProjection:advance(time, inFlds, outFlds)
   local extField = inFlds[1]
   local distf    = outFlds[1]

   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      self.project:advance(t, {}, {distf})
   end
   
   local jacobGeo = extField.geo.jacobGeo
   if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
end

----------------------------------------------------------------------
-- Vlasov-specific VlasovProjection.MaxwellianProjection extends MaxwellianProjection base class.
local MaxwellianProjection = Proto(MaxwellianProjectionParent)

function MaxwellianProjection:lagrangeFix(distf)
   local M0, dM0 = self.species:allocMoment(), self.species:allocMoment()
   local M1      = self.species:allocVectorMoment(self.vdim)
   local dM1     = self.species:allocVectorMoment(self.vdim)
   local M2, dM2 = self.species:allocMoment(), self.species:allocMoment()

   local project = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = function(t,xn) return 0. end,
      onGhosts = true,
   }

   self.species.numDensityCalc:advance(0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   project:setFunc(function(t,xn) return func(t,xn) end)
   project:advance(0.0, {}, {dM0})
   dM0:accumulate(-1.0, M0)

   self.species.momDensityCalc:advance(0.0, {distf}, {M1})
   func = function (t, zn)
      local drifts = self.driftSpeed(t, zn, self.species)
      if self.vdim == 1 then
	 return self.density(t, zn, self.species) * drifts[1]
      elseif self.vdim == 2 then
	 return self.density(t, zn, self.species) * drifts[1], self.density(t, zn, self.species) * drifts[2]
      else
	 return self.density(t, zn, self.species) * drifts[1], self.density(t, zn, self.species) * drifts[2], self.density(t, zn, self.species) * drifts[3]
      end
   end
   project:setFunc(function(t,xn) return func(t,xn) end)
   project:advance(0.0, {}, {dM1})
   dM1:accumulate(-1.0, M1)

   self.species.ptclEnergyCalc:advance(0.0, {distf}, {M2})
   func = function (t, zn)
      local drifts = self.driftSpeed(t, zn, self.species)
      local out = 0.0
      for i = 1, self.vdim do
	 out = out + drifts[i] * drifts[i]
      end
      out = self.density(t, zn, self.species) *
	 (out + self.vdim*self.temperature(t, zn, self.species)/self.species.mass )
      return out
   end
   project:setFunc(function(t,xn) return func(t,xn) end)
   project:advance(0.0, {}, {dM2})
   dM2:accumulate(-1.0, M2)
   
   local lagFix = Updater.LagrangeFix {
      onGrid     = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
      mode       = 'vlasov'
   }
   lagFix:advance(0.0, {dM0, dM1, dM2}, {distf})
end

function MaxwellianProjection:advance(t, inFlds, outFlds)
   local distf = outFlds[1]
   local extField = inFlds[1]
   self.project:advance(t, {}, {distf})
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end
   if self.exactLagFixM012 then
      self:lagrangeFix(distf)
   end

   local jacobGeo = extField.geo.jacobGeo
   if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
end

----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
