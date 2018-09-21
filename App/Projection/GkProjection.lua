-- Gkyl ------------------------------------------------------------------------
--
-- App support code: GkProjection object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

--local Time = require "Lib.Time"
local MaxwellianProjectionParrent = require ("App.Projection.KineticProjection").MaxwellianProjection
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"


local MaxwellianProjection = Proto(MaxwellianProjectionParrent)

function MaxwellianProjection:lagrangeFix(distf)
   local M0, dM0 = self.species:allocMoment(), self.species:allocMoment()
   local M1, dM1 = self.species:allocMoment(), self.species:allocMoment()
   local M2, dM2 = self.species:allocMoment(), self.species:allocMoment()

   self.species.numDensityCalc:advance(0.0, 0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   local project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {dM0})
   dM0:accumulate(-1.0, M0)

   self.species.momDensityCalc:advance(0.0, 0.0, {distf}, {M1})
   func = function (t, zn)
      return self.density(t, zn, self.species) *
	 self.driftSpeed(t, zn, self.species)
   end
   project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {dM1})
   dM1:accumulate(-1.0, M1)


   self.species.ptclEnergyCalc:advance(0.0, 0.0, {distf}, {M2})
   func = function (t, zn)
      local drifts = self.driftSpeed(t, zn, self.species)
      if self.numVelDims == 1 then
	 return self.density(t, zn, self.species) *
	    (self.driftSpeed(t, zn, self.species)*self.driftSpeed(t, zn, self.species) + self.temperature(t, zn, self.species)/self.species.mass )
      else
	 return self.density(t, zn, self.species) *
	    (self.driftSpeed(t, zn, self.species)*self.driftSpeed(t, zn, self.species) + 3*self.temperature(t, zn, self.species)/self.species.mass )
      end
   end
   project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {dM2})
   dM2:accumulate(-1.0, M2)

   local lagFix = Updater.LagrangeFix {
      onGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      mode = 'Gk',
      mass = self.species.mass,
   }
   lagFix:advance(0.0, 0.0, {dM0, dM1, dM2, self.species.bmag}, {distf})
end

function MaxwellianProjection:run(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end
   if self.exactLagFixM012 then
      self:lagrangeFix(distf)
   end
end


----------------------------------------------------------------------
return {
   MaxwellianProjection = MaxwellianProjection,
}
