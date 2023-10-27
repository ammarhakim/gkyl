-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for GkSpecies.
--
-- We assume that integrated diagnostics can depend on grid diagnostics, but
-- not the other way around. We also assume that the integrated diagnostics are
-- always computed before the field diagnostics, since they may be computed more
-- frequently. This allows us to reset the state in calcIntegratedDiagostics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto         = require "Lib.Proto"
local DiagsImplBase = require "App.Diagnostics.DiagnosticsImplBase"
local DataStruct    = require "DataStruct"
local Updater       = require "Updater"
local Constants     = require "Lib.Constants"

local implementation = function()
   -- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
   local _M0 = Proto(DiagsImplBase)
   function _M0:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.fieldwJacobGeo   = owner:allocMoment()  -- So intM0 can depend on M0.
      self.updater          = owner.numDensityCalc or specIn.numDensityCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.fieldwJacobGeo})
      self.divideByJacobGeo(tm, self.fieldwJacobGeo, self.field)
   end

   -- ~~~~ Momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _M1 = Proto(DiagsImplBase)
   function _M1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.fieldwJacobGeo   = owner:allocMoment()  -- So intM1 can depend on M1.
      self.updater          = owner.momDensityCalc or specIn.momDensityCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M1:getType() return "grid" end
   function _M1:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.fieldwJacobGeo})
      self.divideByJacobGeo(tm, self.fieldwJacobGeo, self.field)
   end

   -- ~~~~ Particle energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2 = Proto(DiagsImplBase)
   function _M2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.fieldwJacobGeo   = owner:allocMoment()  -- So intM2 can depend on M2.
      self.updater          = owner.ptclEnergyCalc or specIn.ptclEnergyCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M2:getType() return "grid" end
   function _M2:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.fieldwJacobGeo})
      self.divideByJacobGeo(tm, self.fieldwJacobGeo, self.field)
   end

   -- ~~~~ Mean flow velocity ~~~~~~~~~~~~~~~~~~~~~~
   local _Upar = Proto(DiagsImplBase)
   function _Upar:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakDivide or specIn.confWeakDivide
      self.done    = false
   end
   function _Upar:getType() return "grid" end
   function _Upar:getDependencies() return {"M0","M1"} end
   function _Upar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M1        = diags["M0"].field, diags["M1"].field
      self.updater:advance(tm, {M0, M1}, {self.field})
   end

   -- ~~~~ Mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2Flow = Proto(DiagsImplBase)
   function _M2Flow:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakMultiply or specIn.confWeakMultiply
      self.done    = false
   end
   function _M2Flow:getType() return "grid" end
   function _M2Flow:getDependencies() return {"M1","Upar"} end
   function _M2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1, Upar      = diags["M1"].field, diags["Upar"].field
      self.updater:advance(tm, {M1, Upar}, {self.field})
   end

   -- ~~~~ Thermal energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2Thermal = Proto(DiagsImplBase)
   function _M2Thermal:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _M2Thermal:getType() return "grid" end
   function _M2Thermal:getDependencies() return {"M2","M2Flow"} end
   function _M2Thermal:advance(tm, inFlds, outFlds)
      local diags      = inFlds[2]
      local M2, M2Flow = diags["M2"].field, diags["M2Flow"].field
      self.field:combine(1., M2, -1., M2Flow)
   end

   -- ~~~~ Parallel energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2par = Proto(DiagsImplBase)
   function _M2par:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.updater          = owner.M2parCalc or specIn.M2parCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M2par:getType() return "grid" end
   function _M2par:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
      self.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ mu*B moment, perpendicular energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2perp = Proto(DiagsImplBase)
   function _M2perp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.updater          = owner.M2perpCalc or specIn.M2perpCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M2perp:getType() return "grid" end
   function _M2perp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
      self.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Parallel heat flux density ~~~~~~~~~~~~~~~~~~~~~~
   local _M3par = Proto(DiagsImplBase)
   function _M3par:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.updater          = owner.M3parCalc or specIn.M3parCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M3par:getType() return "grid" end
   function _M3par:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
      self.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Parallel flux of perpendicular energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M3perp = Proto(DiagsImplBase)
   function _M3perp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field            = owner:allocMoment()
      self.updater          = owner.M3perpCalc or specIn.M3perpCalc
      self.divideByJacobGeo = owner.divideByJacobGeo or specIn.divideByJacobGeo
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _M3perp:getType() return "grid" end
   function _M3perp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
      self.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Temp = Proto(DiagsImplBase)
   function _Temp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakDivide or specIn.confWeakDivide
      self.done    = false
   end
   function _Temp:getType() return "grid" end
   function _Temp:getDependencies() return {"M2Thermal","M0"} end
   function _Temp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2Thermal = diags["M0"].field, diags["M2Thermal"].field
      self.updater:advance(tm, {M0, M2Thermal}, {self.field})
      self.field:scale(specIn.mass/specIn.vDegFreedom)
   end

   -- ~~~~ Thermal speed squared ~~~~~~~~~~~~~~~~~~~~~~
   local _VtSq = Proto(DiagsImplBase)
   function _VtSq:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakDivide or specIn.confWeakDivide
      self.done    = false
   end
   function _VtSq:getType() return "grid" end
   function _VtSq:getDependencies() return {"M2Thermal","M0"} end
   function _VtSq:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2Thermal = diags["M0"].field, diags["M2Thermal"].field
      self.updater:advance(tm, {M0, M2Thermal}, {self.field})
      self.field:scale(1./specIn.vDegFreedom)
   end

   -- ~~~~ Parallel temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tpar = Proto(DiagsImplBase)
   function _Tpar:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakDivide or specIn.confWeakDivide
      self.done    = false
   end
   function _Tpar:getType() return "grid" end
   function _Tpar:getDependencies() return {"M2Flow","M2par","M0"} end
   function _Tpar:advance(tm, inFlds, outFlds)
      local specIn, diags     = inFlds[1], inFlds[2]
      local M2Flow, M2par, M0 = diags["M2Flow"].field, diags["M2par"].field, diags["M0"].field
      self.field:combine(1., M2par, -1., M2Flow)
      self.updater:advance(tm, {M0, self.field}, {self.field})
      self.field:scale(specIn.mass)
   end

   -- ~~~~ Perpendicular temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tperp = Proto(DiagsImplBase)
   function _Tperp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakDivide or specIn.confWeakDivide
      self.done    = false
   end
   function _Tperp:getType() return "grid" end
   function _Tperp:getDependencies() return {"M2perp","M0"} end
   function _Tperp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2perp    = diags["M0"].field, diags["M2perp"].field
      self.updater:advance(tm, {M0, M2perp}, {self.field})
      self.field:scale(specIn.mass)
   end

   -- ~~~~ Plasma beta ~~~~~~~~~~~~~~~~~~~~~~
   local _Beta = Proto(DiagsImplBase)
   function _Beta:fullInit(diagApp, specIn, fieldIn, owner)
      self.field     = owner:allocMoment()
      self.bmagInvSq = owner.bmagInvSq
      self.updater   = owner.confWeakMultiply or specIn.confWeakMultiply
      self.done      = false
   end
   function _Beta:getType() return "grid" end
   function _Beta:getDependencies() return {"Temp","M0"} end
   function _Beta:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, Temp      = diags["M0"].field, diags["Temp"].field
      self.updater:advance(tm, {M0, Temp}, {self.field})
      self.updater:advance(tm, {self.bmagInvSq, self.field}, {self.field})
      self.field:scale(2*Constants.MU0)
   end

   -- ~~~~ Hamiltonian energy: total particle energy including potential energy ~~~~~~~~~~~~~~~~~~~~~~
   local _Energy = Proto(DiagsImplBase)
   function _Energy:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocMoment()
      self.updater = owner.confWeakMultiply or specIn.confWeakMultiply
      self.evalPhi = owner.evalOnConfBoundary
         and function(phiIn) return owner:evalOnConfBoundary(phiIn, owner.confBoundaryField) end
         or function(phiIn) return phiIn end
      self.done    = false
   end
   function _Energy:getType() return "grid" end
   function _Energy:getDependencies() return {"M0","M2"} end
   function _Energy:advance(tm, inFlds, outFlds)
      local specIn, diags, field = inFlds[1], inFlds[2], inFlds[3]
      local M0, M2 = diags["M0"].field, diags["M2"].field
      local phi    = self.evalPhi(field:rkStepperFields()[1].phi)
      self.updater:advance(tm, {M0, phi}, {self.field})
      self.field:scale(specIn.charge)
      self.field:accumulate(0.5*specIn.mass, M2)
   end

   -- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intM0:getDependencies() return {"M0"} end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0 = diags["M0"].fieldwJacobGeo
      self.updater:advance(tm, {M0}, {self.field})
   end

   -- ~~~~ Integrated momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM1 = Proto(DiagsImplBase)
   function _intM1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intM1:getDependencies() return {"M1"} end
   function _intM1:getType() return "integrated" end
   function _intM1:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1 = diags["M1"].fieldwJacobGeo
      self.updater:advance(tm, {M1}, {self.field})
   end

   -- ~~~~ Integrated particle kinetic energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2 = Proto(DiagsImplBase)
   function _intM2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intM2:getDependencies() return {"M2"} end
   function _intM2:getType() return "integrated" end
   function _intM2:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["M2"].fieldwJacobGeo
      self.updater:advance(tm, {M2}, {self.field})
   end

   -- ~~~~ Integrated particle kinetic energy density (with mass/2 factor) ~~~~~~~~~~~~~~~~~~~~~~
   local _intKE = Proto(DiagsImplBase)
   function _intKE:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intKE:getDependencies() return {"M2"} end
   function _intKE:getType() return "integrated" end
   function _intKE:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["M2"].fieldwJacobGeo
      self.updater:advance(tm, {M2, 0.5*specIn.mass}, {self.field})
   end

   -- ~~~~ Integrated particle energy density (including potential) ~~~~~~~~~~~~~~~~~~~~~~
   local _intEnergy = Proto(DiagsImplBase)
   function _intEnergy:fullInit(diagApp, specIn, fieldIn, owner)
      self.field              = owner:allocIntMoment()
      self.fieldAux           = owner:allocMoment()
      self.updater            = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.multiplyByJacobGeo = owner.multiplyByJacobGeo and owner.multiplyByJacobGeo or specIn.multiplyByJacobGeo
      self.done = false
   end
   function _intEnergy:getDependencies() return {"Energy"} end
   function _intEnergy:getType() return "integrated" end
   function _intEnergy:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local ptclEnergy    = diags["Energy"].field
      self.multiplyByJacobGeo(tm, ptclEnergy, self.fieldAux)
      self.updater:advance(tm, {self.fieldAux}, {self.field})
   end

   -- ~~~~ Integrated mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Flow = Proto(DiagsImplBase)
   function _intM2Flow:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intM2Flow:getDependencies() return {"M2Flow"} end
   function _intM2Flow:getType() return "integrated" end
   function _intM2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Flow = diags["M2Flow"].field
      self.updater:advance(tm, {M2Flow}, {self.field})
   end

   -- ~~~~ Integrated thermal energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Thermal = Proto(DiagsImplBase)
   function _intM2Thermal:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = owner.volIntegral and owner.volIntegral.scalar or specIn.volIntegral.scalar
      self.done    = false
   end
   function _intM2Thermal:getDependencies() return {"M2Thermal"} end
   function _intM2Thermal:getType() return "integrated" end
   function _intM2Thermal:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Thermal     = diags["M2Thermal"].field
      self.updater:advance(tm, {M2Thermal}, {self.field})
   end

   -- ~~~~ L1 norm (absolute value) of the distribution function ~~~~~~~~~~~~~~~~~~~~~~
   local _intL1 = Proto(DiagsImplBase)
   function _intL1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = Updater.CartFieldIntegratedQuantCalc {
         onGrid = specIn.grid,   numComponents = 1,
         basis  = specIn.basis,  operator      = "abs",
      }
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _intL1:getType() return "integrated" end
   function _intL1:advance(tm, inFlds, outFlds)
      local fIn = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
   end

   -- ~~~~ L2 norm of the distribution function ~~~~~~~~~~~~~~~~~~~~~~
   local _intL2 = Proto(DiagsImplBase)
   function _intL2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field   = owner:allocIntMoment()
      self.updater = Updater.CartFieldIntegratedQuantCalc {
         onGrid = specIn.grid,   numComponents = 1,
         basis  = specIn.basis,  operator      = "sq",
      }
      self:setGetF(specIn.perturbedDiagnostics, owner)  -- self.getF differentiates between f and delta-f.
      self.done = false
   end
   function _intL2:getType() return "integrated" end
   function _intL2:advance(tm, inFlds, outFlds)
      local fIn = self.getF()
      self.updater:advance(tm, {fIn}, {self.field})
   end

   return {
      M0        = _M0,
      M1        = _M1,
      M2        = _M2,
      Upar      = _Upar,
      M2Flow    = _M2Flow,
      M2Thermal = _M2Thermal,
      M2par     = _M2par,
      M2perp    = _M2perp,
      M3par     = _M3par,
      M3perp    = _M3perp,
      Temp      = _Temp,
      VtSq      = _VtSq,
      Tpar      = _Tpar,
      Tperp     = _Tperp,
      Beta      = _Beta,
      Energy       = _Energy,
      intM0        = _intM0,
      intM1        = _intM1,
      intM2        = _intM2,
      intM2Flow    = _intM2Flow,
      intM2Thermal = _intM2Thermal,
      intKE        = _intKE,
      intEnergy    = _intEnergy,
      intL1        = _intL1,
      intL2        = _intL2,
   }
end

return implementation

