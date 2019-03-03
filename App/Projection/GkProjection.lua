-- yl ------------------------------------------------------------------------
--
-- App support code: GkProjection object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

--local Time = require "Lib.Time"
local FunctionProjectionParent = require ("App.Projection.KineticProjection").FunctionProjection
local MaxwellianProjectionParent = require ("App.Projection.KineticProjection").MaxwellianProjection
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"

----------------------------------------------------------------------
-- Gk-specific GkProjection.FunctionProjection includes jacobian factors in initFunc
local FunctionProjection = Proto(FunctionProjectionParent)
function FunctionProjection:run(tProj, distf)
   if self.species.jacobPhaseFunc and self.vdim > 1 then
      local initFuncWithoutJacobian = self.initFunc
      self.initFunc = function (t, xn)
         local xconf = {}
         for d = 1, self.cdim do
            xconf[d] = xn[d]
         end
         local J = self.species.jacobPhaseFunc(t,xconf)
         local f = initFuncWithoutJacobian(t,xn)
         return J*f
      end
   end
   if self.species.jacobGeoFunc then
      local initFuncWithoutJacobian = self.initFunc
      self.initFunc = function (t, xn)
         local xconf = {}
         for d = 1, self.cdim do
            xconf[d] = xn[d]
         end
         local J = self.species.jacobGeoFunc(t,xconf)
         local f = initFuncWithoutJacobian(t,xn)
         return J*f
      end
   end
   -- note: don't use self.project as this does not have jacobian factors in initFunc
   local project = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = self.initFunc,
      projectOnGhosts = true
   }
   project:advance(tProj, {}, {distf})
end

----------------------------------------------------------------------
-- Gk-specific GkProjection.MaxwellianProjection extends MaxwellianProjection base class, including 
-- adding jacobian factors in initFunc
local MaxwellianProjection = Proto(MaxwellianProjectionParent)

function MaxwellianProjection:lagrangeFix(distf)
   local M0, dM0 = self.species:allocMoment(), self.species:allocMoment()
   local M1, dM1 = self.species:allocMoment(), self.species:allocMoment()
   local M2, dM2 = self.species:allocMoment(), self.species:allocMoment()

   self.species.numDensityCalc:advance(0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   local project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, {}, {dM0})
   dM0:accumulate(-1.0, M0)

   self.species.momDensityCalc:advance(0.0, {distf}, {M1})
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
   project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, {}, {dM2})
   dM2:accumulate(-1.0, M2)

   local lagFix = Updater.LagrangeFix {
      onGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      mode = 'Gk',
      mass = self.species.mass,
   }
   lagFix:advance(0.0, {dM0, dM1, dM2, self.species.bmag}, {distf})
end

function MaxwellianProjection:scaleM012(distf)
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
      return vpar^2/2*self.initFunc(t,zn)
   end
   local project2par = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = distf2parFunc,
      projectOnGhosts = true
   }
   project2par:advance(0.0, {}, {distf2par})
   if self.vdim > 1 then 
      local distf2perpFunc = function (t, zn)
         local mu = zn[self.cdim+2]
         return mu*sp.bmagFunc(t,zn)/sp.mass*self.initFunc(t,zn)
      end
      local project2perp = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,
         basis = self.phaseBasis,
         evaluate = distf2perpFunc,
         projectOnGhosts = true
      }
      project2perp:advance(0.0, {}, {distf2perp})
   end

   -- calculate (inexact) moments of initial distribution function
   sp.numDensityCalc:advance(0.0, {distf}, {M0})
   sp.M2parCalc:advance(0.0, {distf}, {M2par})
   if self.vdim > 1 then sp.M2perpCalc:advance(0.0, {distf}, {M2perp}) end

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
   projectM0:advance(0.0, {}, {M0_e})

   local M2func = function (t, zn)
      return self.density(t, zn, sp)*self.temperature(t, zn, sp)/sp.mass
   end
   local projectM2 = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = M2func,
      projectOnGhosts = true,
   }
   projectM2:advance(0.0, {}, {M2_e})

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
   weakDivision:advance(0.0, {M0, M0_e}, {M0_mod})
   -- calculate M2par_mod = M2_e / M2par
   weakDivision:advance(0.0, {M2par, M2_e}, {M2par_mod})
   -- calculate M2perp_mod = M2_e / M2perp
   if self.vdim > 1 then weakDivision:advance(0.0, {M2perp, M2_e}, {M2perp_mod}) end
   -- calculate M0par_mod = M0_e / M2par
   weakDivision:advance(0.0, {M2par, M0_e}, {M0par_mod})
   -- calculate M0par_mod2 = M0 / M2par
   weakDivision:advance(0.0, {M2par, M0}, {M0par_mod2})
   if self.vdim > 1 then 
      -- calculate M0perp_mod = M0_e / M2perp
      weakDivision:advance(0.0, {M2perp, M0_e}, {M0perp_mod})
      -- calculate M0perp_mod2 = M0 / M2perp
      weakDivision:advance(0.0, {M2perp, M0}, {M0perp_mod2})
   end
   
   -- calculate M02par_mod = M0par_mod2 * M2par_mod = (M0/M2par)*(M2_e/M2par)
   weakMultiplicationConf:advance(0.0, {M0par_mod2, M2par_mod}, {M02par_mod})
   -- calculate M02perp_mod = M0perp_mod2 * M2perp_mod = (M0/M2perp)*(M2perp_e/M2perp)
   if self.vdim > 1 then weakMultiplicationConf:advance(0.0, {M0perp_mod2, M2perp_mod}, {M02perp_mod}) end

   -- calculate distf modifiers from combinations of moment modifiers
   if self.vdim==1 then 
      distf0_mod:combine(3/2, M0_mod, -1/2, M2par_mod)
   else 
      distf0_mod:combine(5/2, M0_mod, -1/2, M2par_mod, -1, M2perp_mod)
   end
   distf2par_mod:combine(1, M02par_mod, -1, M0par_mod)
   if self.vdim > 1 then distf2perp_mod:combine(1, M02perp_mod, -1, M0perp_mod) end

   -- calculate distf0 = distf0_mod * distf0
   weakMultiplicationPhase:advance(0.0, {distf0_mod, distf0}, {distf0})
   -- calculate distf2par = distf2par_mod * distf2par
   weakMultiplicationPhase:advance(0.0, {distf2par_mod, distf2par}, {distf2par})
   -- calculate distf2perp = distf2perp_mod * distf2perp
   if self.vdim > 1 then weakMultiplicationPhase:advance(0.0, {distf2perp_mod, distf2perp}, {distf2perp}) end

   -- combine and finish
   distf:combine(1, distf0, 1, distf2par)
   if self.vdim > 1 then distf:accumulate(1, distf2perp) end
end

-- this implementation also corrects Upar, still being tested.
--function MaxwellianProjection:scaleM012(distf)
--   local sp = self.species
--
--   -- initialize maxwellian distribution distf0 = FM, along with 
--   -- distf_vpar = vpar*FM, distf_vpar2 = vpar^2*FM, and distf_muB = mu*B*FM
--   local distf0, distf_vpar, distf_vpar2, distf_muB = sp:allocDistf(), sp:allocDistf(), sp:allocDistf(), sp:allocDistf()
--   distf0:copy(distf)
--   local distf_vparFunc = function (t, zn)
--      local vpar = zn[self.cdim+1]
--      return vpar*self.initFunc(t,zn)
--   end
--   local project_vpar = Updater.ProjectOnBasis {
--      onGrid = self.phaseGrid,
--      basis = self.phaseBasis,
--      evaluate = distf_vparFunc,
--      projectOnGhosts = true
--   }
--   project_vpar:advance(0.0, {}, {distf_vpar})
--   local distf_vpar2Func = function (t, zn)
--      local vpar = zn[self.cdim+1]
--      return vpar^2*self.initFunc(t,zn)
--   end
--   local project_vpar2 = Updater.ProjectOnBasis {
--      onGrid = self.phaseGrid,
--      basis = self.phaseBasis,
--      evaluate = distf_vpar2Func,
--      projectOnGhosts = true
--   }
--   project_vpar2:advance(0.0, {}, {distf_vpar2})
--   if self.vdim == 2 then 
--      local distf_muBFunc = function (t, zn)
--         local mu = zn[self.cdim+2]
--         return mu*sp.bmagFunc(t,zn)*self.initFunc(t,zn)
--      end
--      local project_muB = Updater.ProjectOnBasis {
--         onGrid = self.phaseGrid,
--         basis = self.phaseBasis,
--         evaluate = distf_muBFunc,
--         projectOnGhosts = true
--      }
--      project_muB:advance(0.0, {}, {distf_muB})
--   end
--
--   -- initialize weak multiplication/division operators
--   local weakDivision = Updater.CartFieldBinOp {
--      onGrid = self.confGrid,
--      weakBasis = self.confBasis,
--      operation = "Divide",
--      onGhosts = true,
--   }
--   local weakMultiplicationConf = Updater.CartFieldBinOp {
--      onGrid = self.confGrid,
--      weakBasis = self.confBasis,
--      operation = "Multiply",
--      onGhosts = true,
--   }
--   local weakMultiplicationPhase = Updater.CartFieldBinOp {
--      onGrid = self.phaseGrid,
--      weakBasis = self.phaseBasis,
--      fieldBasis = self.confBasis,
--      operation = "Multiply",
--      onGhosts = true,
--   }
--
--   -- calculate (inexact) moments of initial distribution function
--   local Dens, M1, M2par, M2perp = sp:allocMoment(), sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   sp.numDensityCalc:advance(0.0, {distf}, {Dens})
--   sp.momDensityCalc:advance(0.0, {distf}, {M1})
--   sp.M2parCalc:advance(0.0, {distf}, {M2par})
--   if self.vdim == 2 then sp.M2perpCalc:advance(0.0, {distf}, {M2perp}) end
--
--   -- calculate weak moments
--   local Upar, Upar_sq, Tpar, Tperp = sp:allocMoment(), sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   -- Upar = M1/Dens
--   weakDivision:advance(0.0, {Dens, M1}, {Upar})
--   -- Upar^2
--   weakMultiplicationConf:advance(0.0, {Upar, Upar}, {Upar_sq})
--   -- Tpar = m*(M2par/Dens - M1^2/Dens^2) = m*(M2par/Dens - Upar^2)
--   weakDivision:advance(0.0, {Dens, M2par}, {Tpar})
--   Tpar:accumulate(-1.0, Upar_sq)
--   Tpar:scale(sp.mass)
--   -- Tperp = m*M2perp/Dens
--   if self.vdim == 2 then
--      weakDivision:advance(0.0, {Dens, M2perp}, {Tperp})
--      Tperp:scale(sp.mass)
--   end
--
--   -- initialize exact moments
--   local Dens_e, Upar_e, Temp_e = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   local projectDens = Updater.ProjectOnBasis {
--      onGrid = self.confGrid,
--      basis = self.confBasis,
--      evaluate = function(t, zn) return self.density(t, zn, sp) end,
--      projectOnGhosts = true,
--   }
--   projectDens:advance(0.0, {}, {Dens_e})
--
--   local projectUpar = Updater.ProjectOnBasis {
--      onGrid = self.confGrid,
--      basis = self.confBasis,
--      evaluate = function(t, zn) return self.driftSpeed(t, zn, sp) end,
--      projectOnGhosts = true,
--   }
--   projectUpar:advance(0.0, {}, {Upar_e})
--
--   local projectTemp = Updater.ProjectOnBasis {
--      onGrid = self.confGrid,
--      basis = self.confBasis,
--      evaluate = function(t, zn) return self.temperature(t, zn, sp) end,
--      projectOnGhosts = true,
--   }
--   projectTemp:advance(0.0, {}, {Temp_e})
--
--   local unitField, TparInv, TperpInv = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   local projectUnity = Updater.ProjectOnBasis {
--      onGrid = self.confGrid,
--      basis = self.confBasis,
--      evaluate = function(t, zn) return 1.0 end,
--      projectOnGhosts = true,
--   }
--   projectUnity:advance(0.0, {}, {unitField})
--
--   -- calculate TparInv = 1/Tpar
--   weakDivision:advance(0.0, {Tpar, unitField}, {TparInv})
--   -- calculate TperpInv = 1/Tperp
--   if self.vdim == 2 then weakDivision:advance(0.0, {Tperp, unitField}, {TperpInv}) end
--
--   -- calculate modifier terms for correcting Upar, Tpar, and Tperp:
--   local Dens_mod, Upar_mod, Tpar_mod, Tperp_mod = sp:allocMoment(), sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   local Ppar_mod, Pperp_mod = sp:allocMoment(), sp:allocMoment()
--   local UparDiff, UparDiff_sq = sp:allocMoment(), sp:allocMoment()
--   -- calculate Dens_mod = Dens_e / Dens
--   weakDivision:advance(0.0, {Dens, Dens_e}, {Dens_mod})
--   -- calculate Tpar_mod = Temp_e / Tpar
--   weakDivision:advance(0.0, {Tpar, Temp_e}, {Tpar_mod})
--   -- calculate Tperp_mod = Temp_e / Tperp
--   if self.vdim == 2 then weakDivision:advance(0.0, {Tperp, Temp_e}, {Tperp_mod}) end
--
--   -- calculate Ppar_mod = Dens_e*Temp_e/(Dens*Tpar) = Dens_mod * Tpar_mod
--   weakMultiplicationConf:advance(0.0, {Dens_mod, Tpar_mod}, {Ppar_mod})
--   -- calculate Pperp_mod = Dens_e*Temp_e/(Dens*Tperp) = Dens_mod * Tperp_mod
--   if self.vdim == 2 then weakMultiplicationConf:advance(0.0, {Dens_mod, Tperp_mod}, {Pperp_mod}) end
--
--   -- calculate UparDiff = Upar_e - Upar
--   UparDiff:combine(1.0, Upar_e, -1.0, Upar)
--   -- calculate UparDiff_sq = UparDiff^2 = (Upar_e - Upar)^2
--   weakMultiplicationConf:advance(0.0, {UparDiff, UparDiff}, {UparDiff_sq})
--
--   -- calculate Upar_mod = m*Dens_e/Dens*(Upar_e-Upar)/Tpar
--   weakDivision:advance(0.0, {Tpar, UparDiff}, {Upar_mod}) -- (Upar_e-Upar)/Tpar
--   weakMultiplicationConf:advance(0.0, {Dens_mod, Upar_mod}, {Upar_mod}) -- Dens_e/Dens*(Upar_e-Upar)/Tpar
--   Upar_mod:scale(sp.mass)
--
--   -- calculate Tpar_mod = Dens_e*Tpar_e/(Dens*Tpar) - Dens_e/Dens + m*Dens_e/Dens*(Upar_e-Upar)^2/Tpar
--   weakDivision:advance(0.0, {Tpar, UparDiff_sq}, {Tpar_mod}) -- (Upar_e-Upar)^2/Tpar
--   weakMultiplicationConf:advance(0.0, {Dens_mod, Tpar_mod}, {Tpar_mod}) -- Dens_e/Dens*(Upar_e-Upar)^2/Tpar
--   Tpar_mod:combine(1.0, Ppar_mod, -1.0, Dens_mod, sp.mass, Tpar_mod)
--   
--   -- calculate Tperp_mod = Dens_e/Dens*Tperp_e/Tperp - Dens_e/Dens
--   if self.vdim == 2 then Tperp_mod:combine(1.0, Pperp_mod, -1.0, Dens_mod) end
--
--   -- calculate weighted maxwellians to correct Upar, Tpar, and Tperp:
--   local distf_Upar, distf_Tpar, distf_Tperp = sp:allocDistf(), sp:allocDistf(), sp:allocDistf()
--   -- calculate distf_Upar = (vpar - Upar)*FM
--   weakMultiplicationPhase:advance(0.0, {Upar, distf0}, {distf_Upar}) -- Upar*FM
--   distf_Upar:combine(1.0, distf_vpar, -1.0, distf_Upar) 
--
--   -- calculate distf_Tpar = (m*(vpar - Upar)^2/(2*Tpar) - 1/2)*FM
--   weakMultiplicationPhase:advance(0.0, {Upar, distf0}, {distf_Tpar}) -- Upar*FM
--   distf_Tpar:combine(1.0, distf_vpar, -.5, distf_Tpar) -- vpar*FM - Upar/2*FM
--   weakMultiplicationPhase:advance(0.0, {Upar, distf_Tpar}, {distf_Tpar}) -- Upar*(vpar*FM - Upar/2*FM)
--   distf_Tpar:combine(.5, distf_vpar2, -1, distf_Tpar) -- vpar^2/2*FM - Upar*(vpar*FM - Upar/2*FM)
--   weakMultiplicationPhase:advance(0.0, {TparInv, distf_Tpar}, {distf_Tpar}) -- (vpar^2/2*FM - Upar*(vpar*FM - Upar/2*FM))/Tpar
--   distf_Tpar:combine(sp.mass, distf_Tpar, -.5, distf0)
--  
--   -- calculate distf_Tperp = (mu*B/Tperp - 1)*FM
--   if self.vdim == 2 then 
--      weakMultiplicationPhase:advance(0.0, {TperpInv, distf_muB}, {distf_Tperp}) -- mu*B/Tperp*FM
--      distf_Tperp:combine(1.0, distf_Tperp, -1.0, distf0)
--   end
--
--   -- multiply weighted maxwellians by modifier terms:
--   -- calculate distf0 = Dens_mod*distf0
--   weakMultiplicationPhase:advance(0.0, {Dens_mod, distf0}, {distf0})
--   -- calculate distf_Upar = Upar_mod*distf_Upar
--   weakMultiplicationPhase:advance(0.0, {Upar_mod, distf_Upar}, {distf_Upar})
--   -- calculate distf_Tpar = Tpar_mod*distf_Tpar
--   weakMultiplicationPhase:advance(0.0, {Tpar_mod, distf_Tpar}, {distf_Tpar})
--   -- calculate distf_Tperp = Tperp_mod*distf_Tperp
--   if self.vdim == 2 then weakMultiplicationPhase:advance(0.0, {Tperp_mod, distf_Tperp}, {distf_Tperp}) end
--
--   -- combine and finish
--   distf:combine(1, distf0, 1, distf_Upar, distf_Tpar)
--   if self.vdim == 2 then distf:accumulate(1, distf_Tperp) end
--end

function MaxwellianProjection:run(tProj, distf)
   if self.species.jacobPhaseFunc and self.vdim > 1 then
      local initFuncWithoutJacobian = self.initFunc
      self.initFunc = function (t, xn)
         local xconf = {}
         for d = 1, self.cdim do
            xconf[d] = xn[d]
         end
         local J = self.species.jacobPhaseFunc(t,xconf)
         local f = initFuncWithoutJacobian(t,xn)
         return J*f
      end
   end
   -- for geometry jacobian, scale density function so that jacobian factor
   -- is retained even after rescaling distf
   if self.species.jacobGeoFunc then
      local densityWithoutJacobian = self.density
      self.density = function (t, xn, sp)
         local xconf = {}
         for d = 1, self.cdim do
            xconf[d] = xn[d]
         end
         local J = self.species.jacobGeoFunc(t,xconf)
         local n = densityWithoutJacobian(t,xn,sp)
         return J*n
      end
   end
   -- note: don't use self.project as this does not have jacobian factors in initFunc
   local project = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = self.initFunc,
      projectOnGhosts = true
   }
   project:advance(tProj, {}, {distf})
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
   FunctionProjection = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
