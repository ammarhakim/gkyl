-- Gkyl ------------------------------------------------------------------------
--
-- App support code: GkProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local lume       = require "Lib.lume"
local Updater    = require "Updater"
local xsys       = require "xsys"
local DataStruct = require "DataStruct"
local FunctionProjectionParent   = require ("App.Projection.KineticProjection").FunctionProjection
local MaxwellianProjectionParent = require ("App.Projection.KineticProjection").MaxwellianProjection
local BiMaxwellianProjectionParent = require ("App.Projection.KineticProjection").BiMaxwellianProjection

--------------------------------------------------------------------------------
-- Gk-specific GkProjection.FunctionProjection includes Jacobian factors in initFunc.
local FunctionProjection = Proto(FunctionProjectionParent)
function FunctionProjection:advance(time, inFlds, outFlds)
   local extField = inFlds[1]
   local distf    = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      if self.species.jacobPhaseFunc and self.vdim > 1 then
         local initFuncWithoutJacobian = self.initFunc
         self.initFunc = function (t, xn)
            local xconf = {}
            for d = 1, self.cdim do xconf[d] = xn[d] end
            local J = self.species.jacobPhaseFunc(t,xconf)
            local f = initFuncWithoutJacobian(t,xn)
            return J*f
         end
      end

      -- Note: don't use self.project as this does not have jacobian factors in initFunc.
      local project = Updater.ProjectOnBasis {
         onGrid   = self.phaseGrid,
         basis    = self.phaseBasis,
         evaluate = self.initFunc,
         onGhosts = true
      }
      project:advance(time, {}, {distf})
   end

   local jacobGeo = extField.geo.jacobGeo
   if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
end

--------------------------------------------------------------------------------
-- Gk-specific GkProjection.MaxwellianProjection extends MaxwellianProjection base class, including 
-- adding jacobian factors in initFunc.
local MaxwellianProjection = Proto(MaxwellianProjectionParent)

function MaxwellianProjection:allocConfField()
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
      onGrid     = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
      mode       = 'gk',
      mass       = self.species.mass,
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
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts  = true,
   }
   local weakMultiplicationConf = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Multiply",
      onGhosts  = true,
   }
   local weakMultiplicationPhase = Updater.CartFieldBinOp {
      onGrid     = self.phaseGrid,
      weakBasis  = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
      onGhosts   = true,
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
--      onGrid   = self.phaseGrid,
--      basis    = self.phaseBasis,
--      evaluate = distf_vparFunc,
--      onGhosts = true
--   }
--   project_vpar:advance(0.0, {}, {distf_vpar})
--   local distf_vpar2Func = function (t, zn)
--      local vpar = zn[self.cdim+1]
--      return vpar^2*self.initFunc(t,zn)
--   end
--   local project_vpar2 = Updater.ProjectOnBasis {
--      onGrid   = self.phaseGrid,
--      basis    = self.phaseBasis,
--      evaluate = distf_vpar2Func,
--      onGhosts = true
--   }
--   project_vpar2:advance(0.0, {}, {distf_vpar2})
--   if self.vdim == 2 then 
--      local distf_muBFunc = function (t, zn)
--         local mu = zn[self.cdim+2]
--         return mu*sp.bmagFunc(t,zn)*self.initFunc(t,zn)
--      end
--      local project_muB = Updater.ProjectOnBasis {
--         onGrid   = self.phaseGrid,
--         basis    = self.phaseBasis,
--         evaluate = distf_muBFunc,
--         onGhosts = true
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
--      onGrid   = self.confGrid,
--      basis    = self.confBasis,
--      evaluate = function(t, zn) return self.density(t, zn, sp) end,
--      onGhosts = true,
--   }
--   projectDens:advance(0.0, {}, {Dens_e})
--
--   local projectUpar = Updater.ProjectOnBasis {
--      onGrid   = self.confGrid,
--      basis    = self.confBasis,
--      evaluate = function(t, zn) return self.driftSpeed(t, zn, sp) end,
--      onGhosts = true,
--   }
--   projectUpar:advance(0.0, {}, {Upar_e})
--
--   local projectTemp = Updater.ProjectOnBasis {
--      onGrid   = self.confGrid,
--      basis    = self.confBasis,
--      evaluate = function(t, zn) return self.temperature(t, zn, sp) end,
--      onGhosts = true,
--   }
--   projectTemp:advance(0.0, {}, {Temp_e})
--
--   local unitField, TparInv, TperpInv = sp:allocMoment(), sp:allocMoment(), sp:allocMoment()
--   local projectUnity = Updater.ProjectOnBasis {
--      onGrid   = self.confGrid,
--      basis    = self.confBasis,
--      evaluate = function(t, zn) return 1.0 end,
--      onGhosts = true,
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

function MaxwellianProjection:advance(time, inFlds, outFlds)
   local extField = inFlds[1]
   local distf    = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      local bmag = extField.geo.bmag
      -- Project the moments onto configuration-space basis.
      local confProject = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 0. end,   -- Set below.
         basis  = self.confBasis,  onGhosts = true
      }
      local numDens = self:allocConfField()
      local uPar    = self:allocConfField()
      local vtSq    = self:allocConfField()
      --Create Open BC Updaters for Lower and Upper Edges
      local nonPeriodicBCs = {}   -- List of non-periodic BCs to apply.
      -- Function to construct a BC updater.
      local function makeOpenBcUpdater(dir, edge)
         local bcOpen = function(dir, tm, idxIn, fIn, fOut)
            -- Requires skinLoop = "pointwise".
            self.confBasis:flipSign(dir, fIn, fOut)
         end
         return Updater.Bc {
            onGrid = self.confGrid,  edge = edge,
            dir    = dir,        boundaryConditions = {bcOpen},
            skinLoop = "pointwise",
         }
      end
      -- For non-periodic dirs, use open BCs. It's the most sensible choice given that the
      -- coordinate mapping could diverge outside of the interior domain.
      for dir = 1, self.cdim do
         if not lume.any(self.confGrid:getPeriodicDirs(), function(t) return t==dir end) then
         nonPeriodicBCs["lower" .. tostring(dir)] = makeOpenBcUpdater(dir, "lower")
         nonPeriodicBCs["upper" .. tostring(dir)] = makeOpenBcUpdater(dir, "upper")
         end
      end
      --Function to Add and Apply Open BCs to a confField
      local function applyOpenBcs(fld)
         for _, bc in pairs(nonPeriodicBCs) do
            bc:advance(tCurr, {}, {fld})
         end
         fld:sync(true)
      end

      if self.densityFromFile then
         self.confFieldIo:read(numDens, self.density,false)
         applyOpenBcs(numDens)
      else
        if self.exactScaleM0 then
           -- Use a unit density because we are going to rescale the density anyways,
           -- and it is easier to weak-divide by something close to unity.
           confProject:setFunc(function(t, xn) return 1. end)
        else
           confProject:setFunc(self.density)
        end
        confProject:advance(time, {}, {numDens})
      end
      if self.driftSpeedFromFile then
         self.confFieldIo:read(uPar, self.driftSpeed,false)
         applyOpenBcs(uPar)
      else
         confProject:setFunc(self.driftSpeed)
         confProject:advance(time, {}, {uPar})
      end
      if self.temperatureFromFile then
         self.confFieldIo:read(vtSq, self.temperature,false)
         applyOpenBcs(vtSq)
      else
         confProject:setFunc(self.temperature)
         confProject:advance(time, {}, {vtSq})
      end
      vtSq:scale(1./self.mass)
      -- Project the Maxwellian. It includes a factor of jacobPhase=B*_||.
      local projMaxwell = Updater.MaxwellianOnBasis {
         onGrid     = self.phaseGrid,   confBasis = self.confBasis,
         phaseBasis = self.phaseBasis,  mass      = self.mass,
         confGrid   = self.confGrid,    onGhosts  = true,
      }
      projMaxwell:advance(time,{numDens,uPar,vtSq,bmag},{distf})
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

--------------------------------------------------------------------------------
-- Gk-specific GkProjection.BiMaxwellianProjection extends BiMaxwellianProjection base class, including 
-- adding jacobian factors in initFunc.
local BiMaxwellianProjection = Proto(BiMaxwellianProjectionParent)

function BiMaxwellianProjection:allocConfField()
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

function BiMaxwellianProjection:advance(time, inFlds, outFlds)
   local extField = inFlds[1]
   local distf    = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      local bmag = extField.geo.bmag
      -- Project the moments onto configuration-space basis.
      local confProject = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 0. end,   -- Set below.
         basis  = self.confBasis,  onGhosts = true
      }
      local numDens  = self:allocConfField()
      local uPar     = self:allocConfField()
      local vtparSq  = self:allocConfField()
      local vtperpSq = self:allocConfField()
      if self.exactScaleM0 then
         -- Use a unit density because we are going to rescale the density anyways,
         -- and it is easier to weak-divide by something close to unity.
         confProject:setFunc(function(t, xn) return 1. end)
      else
         confProject:setFunc(self.density)
      end
      confProject:advance(time, {}, {numDens})
      confProject:setFunc(self.driftSpeed)
      confProject:advance(time, {}, {uPar})
      confProject:setFunc(self.parallelTemperature)
      confProject:advance(time, {}, {vtparSq})
      confProject:setFunc(self.perpendicularTemperature)
      confProject:advance(time, {}, {vtperpSq})
      vtparSq:scale(1./self.mass)
      vtperpSq:scale(1./self.mass)
      -- Project the BiMaxwellian. It includes a factor of jacobPhase=B*_||.
      local projMaxwell = Updater.BiMaxwellianOnBasis {
         onGrid     = self.phaseGrid,   confBasis = self.confBasis,
         phaseBasis = self.phaseBasis,  mass      = self.mass,
         confGrid   = self.confGrid,    onGhosts  = true,
      }
      projMaxwell:advance(time,{numDens,uPar,vtparSq,vtperpSq,bmag},{distf})
   end

   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end

   local jacobGeo = extField.geo.jacobGeo
   if jacobGeo then self.weakMultiplyConfPhase:advance(0, {distf, jacobGeo}, {distf}) end
end


----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
   BiMaxwellianProjection = BiMaxwellianProjection,
}
