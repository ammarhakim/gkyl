-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for VlasovSpecies.
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
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.numDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _M1 = Proto(DiagsImplBase)
   function _M1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M1:getType() return "grid" end
   function _M1:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.momDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Particle energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2 = Proto(DiagsImplBase)
   function _M2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M2:getType() return "grid" end
   function _M2:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.ptclEnergyCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Mean flow velocity ~~~~~~~~~~~~~~~~~~~~~~
   local _uPar = Proto(DiagsImplBase)
   function _uPar:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _uPar:getType() return "grid" end
   function _uPar:getDependencies() return {"M0","M1"} end
   function _uPar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M1        = diags["M0"].field, diags["M1"].field
      specIn.confWeakDivide:advance(tm, {M0, M1}, {self.field})
   end

   -- ~~~~ Mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2Flow = Proto(DiagsImplBase)
   function _M2Flow:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _M2Flow:getType() return "grid" end
   function _M2Flow:getDependencies() return {"M1","uPar"} end
   function _M2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1, uPar      = diags["M1"].field, diags["uPar"].field
      specIn.confWeakMultiply:advance(tm, {M1, uPar}, {self.field})
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
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M2par:getType() return "grid" end
   function _M2par:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.M2parCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ mu*B moment, perpendicular energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2perp = Proto(DiagsImplBase)
   function _M2perp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M2perp:getType() return "grid" end
   function _M2perp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.M2perpCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Parallel heat flux density ~~~~~~~~~~~~~~~~~~~~~~
   local _M3par = Proto(DiagsImplBase)
   function _M3par:fullInit(diagApp, specIn, fieldIn, owner)
      self.field    = owner:allocMoment()
      self.owner    = owner
      self.updaters = Updater.DistFuncMomentCalc {
         onGrid     = specIn.grid,    confBasis = specIn.confBasis,
         phaseBasis = specIn.basis,   moment    = "GkM3par",
         gkfacs     = {specIn.mass, specIn.bmag},
      }
      self.done = false
   end
   function _M3par:getType() return "grid" end
   function _M3par:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Parallel flux of perpendicular energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M3perp = Proto(DiagsImplBase)
   function _M3perp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field    = owner:allocMoment()
      self.owner    = owner
      self.updaters = Updater.DistFuncMomentCalc {
         onGrid     = specIn.grid,    confBasis = specIn.confBasis,
         phaseBasis = specIn.basis,   moment    = "GkM3perp",
         gkfacs     = {specIn.mass, specIn.bmag},
      }
      self.done = false
   end
   function _M3perp:getType() return "grid" end
   function _M3perp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Temp = Proto(DiagsImplBase)
   function _Temp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _Temp:getType() return "grid" end
   function _Temp:getDependencies() return {"M2Thermal","M0"} end
   function _Temp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2Thermal = diags["M0"].field, diags["M2Thermal"].field
      specIn.confWeakDivide:advance(tm, {M0, M2Thermal}, {self.field})
      self.field:scale(specIn.mass/specIn.vDegFreedom)
   end

   -- ~~~~ Thermal speed squared ~~~~~~~~~~~~~~~~~~~~~~
   local _vtSq = Proto(DiagsImplBase)
   function _vtSq:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _vtSq:getType() return "grid" end
   function _vtSq:getDependencies() return {"M2Thermal","M0"} end
   function _vtSq:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2Thermal = diags["M0"].field, diags["M2Thermal"].field
      specIn.confWeakDivide:advance(tm, {M0, M2Thermal}, {self.field})
      self.field:scale(1./specIn.vDegFreedom)
   end

   -- ~~~~ Parallel temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tpar = Proto(DiagsImplBase)
   function _Tpar:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _Tpar:getType() return "grid" end
   function _Tpar:getDependencies() return {"M2Flow","M2par","M0"} end
   function _Tpar:advance(tm, inFlds, outFlds)
      local specIn, diags     = inFlds[1], inFlds[2]
      local M2Flow, M2par, M0 = diags["M2Flow"].field, diags["M2par"].field, diags["M0"].field
      self.field:combine(1., M2par, -1., M2Flow)
      specIn.confWeakDivide:advance(tm, {M0, self.field}, {self.field})
      self.field:scale(specIn.mass)
   end

   -- ~~~~ Perpendicular temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tperp = Proto(DiagsImplBase)
   function _Tperp:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _Tperp:getType() return "grid" end
   function _Tperp:getDependencies() return {"M2perp","M0"} end
   function _Tperp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2perp    = diags["M0"].field, diags["M2perp"].field
      specIn.confWeakDivide:advance(tm, {M0, M2perp}, {self.field})
      self.field:scale(specIn.mass)
   end

   -- ~~~~ plasma beta ~~~~~~~~~~~~~~~~~~~~~~
   local _beta = Proto(DiagsImplBase)
   function _beta:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _beta:getType() return "grid" end
   function _beta:getDependencies() return {"Temp","M0"} end
   function _beta:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, Temp      = diags["M0"].field, diags["Temp"].field
      specIn.confWeakMultiply:advance(tm, {M0, Temp}, {self.field})
      specIn.confWeakMultiply:advance(tm, {specIn.bmagInvSq, self.field}, {self.field})
      self.field:scale(2*Constants.MU0)
   end

   -- ~~~~ Total particle energy (including potential energy) ~~~~~~~~~~~~~~~~~~~~~~
   local _particleEnergy = Proto(DiagsImplBase)
   function _particleEnergy:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = owner:allocMoment()
      self.phi   = fieldIn.potentials[1].phi
      self.done  = false
   end
   function _particleEnergy:getType() return "grid" end
   function _particleEnergy:getDependencies() return {"M0","M2"} end
   function _particleEnergy:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2        = diags["M0"].field, diags["M2"].field
      specIn.confWeakMultiply:advance(tm, {M0, self.phi}, {self.field})
      self.field:scale(specIn.charge)
      self.field:accumulate(0.5*specIn.mass, M2)
   end

   -- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM0:getDependencies() return {"GkM0"} end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0 = diags["GkM0"].field
      specIn.volIntegral:advance(tm, {M0}, {self.field})
   end

   -- ~~~~ Integrated momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM1 = Proto(DiagsImplBase)
   function _intM1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM1:getDependencies() return {"GkM1"} end
   function _intM1:getType() return "integrated" end
   function _intM1:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1i = diags["GkM1"].field
      specIn.volIntegral:advance(tm, {M1i}, {self.field})
   end

   -- ~~~~ Integrated particle kinetic energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2 = Proto(DiagsImplBase)
   function _intM2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM2:getDependencies() return {"GkM2"} end
   function _intM2:getType() return "integrated" end
   function _intM2:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["GkM2"].field
      specIn.volIntegral:advance(tm, {M2}, {self.field})
   end

   -- ~~~~ Integrated particle kinetic energy density (with mass/2 factor) ~~~~~~~~~~~~~~~~~~~~~~
   local _intKE = Proto(DiagsImplBase)
   function _intKE:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intKE:getDependencies() return {"GkM2"} end
   function _intKE:getType() return "integrated" end
   function _intKE:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["GkM2"].field
      specIn.volIntegral:advance(tm, {M2, 0.5*specIn.mass}, {self.field})
   end

   -- ~~~~ Integrated particle energy density (including potential) ~~~~~~~~~~~~~~~~~~~~~~
   local _intHE = Proto(DiagsImplBase)
   function _intHE:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intHE:getDependencies() return {"particleEnergy"} end
   function _intHE:getType() return "integrated" end
   function _intHE:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local ptclEnergy    = diags["particleEnergy"].field
      specIn.volIntegral:advance(tm, {ptclEnergy}, {self.field})
   end

   -- ~~~~ Integrated mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Flow = Proto(DiagsImplBase)
   function _intM2Flow:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM2Flow:getDependencies() return {"M2Flow"} end
   function _intM2Flow:getType() return "integrated" end
   function _intM2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Flow = diags["M2Flow"].field
      specIn.volIntegral:advance(tm, {M2Flow}, {self.field})
   end

   -- ~~~~ Integrated thermal energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Thermal = Proto(DiagsImplBase)
   function _intM2Thermal:fullInit(diagApp, specIn, fieldIn, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM2Thermal:getDependencies() return {"M2Thermal"} end
   function _intM2Thermal:getType() return "integrated" end
   function _intM2Thermal:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Thermal = diags["M2Thermal"].field
      specIn.volIntegral:advance(tm, {M2Thermal}, {self.field})
   end

   -- ~~~~ L1 norm (absolute value) of the distribution function ~~~~~~~~~~~~~~~~~~~~~~
   local _intL1 = Proto(DiagsImplBase)
   function _intL1:fullInit(diagApp, specIn, fieldIn, owner)
      self.field    = DataStruct.DynVector { numComponents = 1 }
      self.updaters = Updater.CartFieldIntegratedQuantCalc {
         onGrid = specIn.grid,   numComponents = 1,
         basis  = specIn.basis,  quantity      = "AbsV",
      }
      self.owner = owner
      self.done  = false
   end
   function _intL1:getType() return "integrated" end
   function _intL1:advance(tm, inFlds, outFlds)
      local fIn = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
   end

   -- ~~~~ L2 norm of the distribution function ~~~~~~~~~~~~~~~~~~~~~~
   local _intL2 = Proto(DiagsImplBase)
   function _intL2:fullInit(diagApp, specIn, fieldIn, owner)
      self.field    = DataStruct.DynVector { numComponents = 1 }
      self.updaters = Updater.CartFieldIntegratedQuantCalc {
         onGrid = specIn.grid,   numComponents = 1,
         basis  = specIn.basis,  quantity      = "V2",
      }
      self.owner = owner
      self.done  = false
   end
   function _intL2:getType() return "integrated" end
   function _intL2:advance(tm, inFlds, outFlds)
      local fIn = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
   end

   return {
      GkM0      = _M0,
      GkM1      = _M1i,
      GkM2      = _M2,
      uPar      = _uPar,
      M2Flow    = _M2Flow,
      M2Thermal = _M2Thermal,
      GkM2par   = _M2par,
      GkM2perp  = _M2perp,
      GkM3par   = _M3par,
      GkM3perp  = _M3par,
      GkTemp    = _Temp,
      vtSq      = _vtSq,
      GkTpar    = _Tpar,
      GkTperp   = _Tperp,
      GkBeta    = _beta,
      GkEnergy  = _particleEnergy,
      intM0        = _intM0,
      intM1        = _intM1,
      intM2        = _intM2,
      intM2Flow    = _intM2Flow,
      intM2Thermal = _intM2Thermal,
      intKE        = _intKE,
      intHE        = _intHE,
      intL1        = _intL1,
      intL2        = _intL2,
   }
end

return implementation

