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
      if self.vdim == 1 then
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

function MaxwellianProjection:scaleM012(distf)
   assert(self.vdim == 2, "scaleM012 currently only implemented for vdim=2")

   local sp = self.species
   local M0, M2par, M2perp = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
   local M0_e, M2_e = sp:allocMoment(), sp:allocMoment()
   local M0_mod, M2par_mod, M2perp_mod = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
   local M0par_mod, M0perp_mod = sp:allocMoment(), sp:allocMoment() 
   local M0par_mod2, M0perp_mod2 = sp:allocMoment(), sp:allocMoment() 
   local M02par_mod, M02perp_mod = sp:allocMoment(), sp:allocMoment() 
   local distf0_mod, distf2par_mod, distf2perp_mod = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()

   -- initialize maxwellian distribution distf0 = FM, along with 
   -- distf2par = m*vpar^2/2*FM and distf2perp = mu*B*FM
   local distf0, distf2par, distf2perp = sp:allocDistf(), sp:allocDistf(), sp:allocDistf()
   distf0:copy(distf)
   local distf2parFunc = function (t, zn)
      local vpar = zn[self.cdim+1]
      return vpar^2/2*sp:Maxwellian(zn, self.density(t, zn, sp),
				self.temperature(t, zn, sp),
				self.driftSpeed(t, zn, sp))
   end
   local project2par = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = distf2parFunc,
      projectOnGhosts = true
   }
   project2par:advance(0.0, 0.0, {}, {distf2par})
   local distf2perpFunc = function (t, zn)
      local mu = zn[self.cdim+2]
      return mu*sp.bmagFunc(t,zn)/sp.mass*sp:Maxwellian(zn, self.density(t, zn, sp),
				self.temperature(t, zn, sp),
				self.driftSpeed(t, zn, sp))
   end
   local project2perp = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = distf2perpFunc,
      projectOnGhosts = true
   }
   project2perp:advance(0.0, 0.0, {}, {distf2perp})

   -- calculate (inexact) moments of initial distribution function
   sp.numDensityCalc:advance(0.0, 0.0, {distf}, {M0})
   sp.M2parCalc:advance(0.0, 0.0, {distf}, {M2par})
   sp.M2perpCalc:advance(0.0, 0.0, {distf}, {M2perp})

   -- initialize exact moments
   local M0func = function (t, zn)
      return self.density(t, zn, sp)
   end
   local projectM0 = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = M0func,
      projectOnGhosts = true,
   }
   projectM0:advance(0.0, 0.0, {}, {M0_e})

   local M2func = function (t, zn)
      return self.density(t, zn, sp)*self.temperature(t, zn, sp)/sp.mass
   end
   local projectM2 = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = M2func,
      projectOnGhosts = true,
   }
   projectM2:advance(0.0, 0.0, {}, {M2_e})

   -- initialize weak multiplication/division operators
   local weakDivision = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts = true,
   }
   local weakMultiplicationConf = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Multiply",
      onGhosts = true,
   }
   local weakMultiplicationPhase = Updater.CartFieldBinOp {
      onGrid = self.phaseGrid,
      weakBasis = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation = "Multiply",
      onGhosts = true,
   }

   -- calculate M0_mod = M0_e / M0
   weakDivision:advance(0.0, 0.0, {M0, M0_e}, {M0_mod})
   -- calculate M2par_mod = M2_e / M2par
   weakDivision:advance(0.0, 0.0, {M2par, M2_e}, {M2par_mod})
   -- calculate M2perp_mod = M2_e / M2perp
   weakDivision:advance(0.0, 0.0, {M2perp, M2_e}, {M2perp_mod})
   -- calculate M0par_mod = M0_e / M2par
   weakDivision:advance(0.0, 0.0, {M2par, M0_e}, {M0par_mod})
   -- calculate M0par_mod2 = M0 / M2par
   weakDivision:advance(0.0, 0.0, {M2par, M0}, {M0par_mod2})
   -- calculate M0perp_mod = M0_e / M2perp
   weakDivision:advance(0.0, 0.0, {M2perp, M0_e}, {M0perp_mod})
   -- calculate M0perp_mod2 = M0 / M2perp
   weakDivision:advance(0.0, 0.0, {M2perp, M0}, {M0perp_mod2})
   
   -- calculate M02par_mod = M0par_mod2 * M2par_mod = (M0/M2par)*(M2_e/M2par)
   weakMultiplicationConf:advance(0.0, 0.0, {M0par_mod2, M2par_mod}, {M02par_mod})
   -- calculate M02perp_mod = M0perp_mod2 * M2perp_mod = (M0/M2perp)*(M2perp_e/M2perp)
   weakMultiplicationConf:advance(0.0, 0.0, {M0perp_mod2, M2perp_mod}, {M02perp_mod})

   -- calculate distf modifiers from combinations of moment modifiers
   distf0_mod:combine(5/2, M0_mod, -1/2, M2par_mod, -1, M2perp_mod)
   distf2par_mod:combine(1, M02par_mod, -1, M0par_mod)
   distf2perp_mod:combine(1, M02perp_mod, -1, M0perp_mod)

   -- calculate distf0 = distf0_mod * distf0
   weakMultiplicationPhase:advance(0.0, 0.0, {distf0_mod, distf0}, {distf0})
   -- calculate distf2par = distf2par_mod * distf2par
   weakMultiplicationPhase:advance(0.0, 0.0, {distf2par_mod, distf2par}, {distf2par})
   -- calculate distf2perp = distf2perp_mod * distf2perp
   weakMultiplicationPhase:advance(0.0, 0.0, {distf2perp_mod, distf2perp}, {distf2perp})

   -- combine and finish
   distf:combine(1, distf0, 1, distf2par, 1, distf2perp)
end

function MaxwellianProjection:run(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   elseif self.exactScaleM012 then
      self:scaleM012(distf)
   end
   if self.exactLagFixM012 then
      self:lagrangeFix(distf)
   end
end


----------------------------------------------------------------------
return {
   MaxwellianProjection = MaxwellianProjection,
}
