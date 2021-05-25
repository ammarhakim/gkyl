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

local implementation = function()
   -- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
   local _M0 = Proto(DiagsImplBase)
   function _M0:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.numDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _M1i = Proto(DiagsImplBase)
   function _M1i:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocVectorMoment(specIn.vdim)
      self.owner = owner
      self.done  = false
   end
   function _M1i:getType() return "grid" end
   function _M1i:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.momDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Particle energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2 = Proto(DiagsImplBase)
   function _M2:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M2:getType() return "grid" end
   function _M2:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.ptclEnergyCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Mean flow velocity ~~~~~~~~~~~~~~~~~~~~~~
   local _Udrift = Proto(DiagsImplBase)
   function _Udrift:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocVectorMoment(specIn.vdim)
      self.done  = false
   end
   function _Udrift:getType() return "grid" end
   function _Udrift:getDependencies() return {"M0","M1i"} end
   function _Udrift:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M1i       = diags["M0"].field, diags["M1i"].field
      specIn.confWeakDivide:advance(tm, {M0, M1i}, {self.field})
   end

   -- ~~~~ Mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2Flow = Proto(DiagsImplBase)
   function _M2Flow:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _M2Flow:getType() return "grid" end
   function _M2Flow:getDependencies() return {"M1i","Udrift"} end
   function _M2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1i, Udrift   = diags["M1i"].field, diags["Udrift"].field
      specIn.confWeakDotProduct:advance(tm, {M1i, Udrift}, {self.field})
   end

   -- ~~~~ Thermal energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2Thermal = Proto(DiagsImplBase)
   function _M2Thermal:fullInit(diagApp, specIn, field, owner)
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

   -- ~~~~ Square velocity tensor density ~~~~~~~~~~~~~~~~~~~~~~
   local _M2ij = Proto(DiagsImplBase)
   function _M2ij:fullInit(diagApp, specIn, field, owner)
      self.field    = owner:allocVectorMoment(specIn.vdim*(specIn.vdim+1)/2)
      self.owner    = owner
      self.updaters = Updater.DistFuncMomentCalc {
         onGrid     = specIn.grid,    confBasis = specIn.confBasis,
         phaseBasis = specIn.basis,   moment    = "M2ij",
      }
      self.done = false
   end
   function _M2ij:getType() return "grid" end
   function _M2ij:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Diagonal heat flux density ~~~~~~~~~~~~~~~~~~~~~~
   local _M3i = Proto(DiagsImplBase)
   function _M3i:fullInit(diagApp, specIn, field, owner)
      self.field    = owner:allocVectorMoment(specIn.vdim)
      self.owner    = owner
      self.updaters = Updater.DistFuncMomentCalc {
         onGrid     = specIn.grid,    confBasis = specIn.confBasis,
         phaseBasis = specIn.basis,   moment    = "M3i",
      }
      self.done = false
   end
   function _M3i:getType() return "grid" end
   function _M3i:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      self.updaters:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field, self.field)
   end

   -- ~~~~ Thermal speed squared ~~~~~~~~~~~~~~~~~~~~~~
   local _VtSq = Proto(DiagsImplBase)
   function _VtSq:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocMoment()
      self.done  = false
   end
   function _VtSq:getType() return "grid" end
   function _VtSq:getDependencies() return {"M2Thermal","M0"} end
   function _VtSq:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0, M2Thermal = diags["M0"].field, diags["M2Thermal"].field
      specIn.confWeakDivide:advance(tm, {M0, M2Thermal}, {self.field})
      self.field:scale(1./specIn.vdim)
   end

   -- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
      self:setVolIntegral(specIn)
   end
   function _intM0:getDependencies() return {"M0"} end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0 = diags["M0"].field
      specIn.volIntegral.comps1:advance(tm, {M0}, {self.field})
   end

   -- ~~~~ Integrated momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM1i = Proto(DiagsImplBase)
   function _intM1i:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = specIn.vdim }
      self.done  = false
      self:setVolIntegral(specIn)
   end
   function _intM1i:getDependencies() return {"M1i"} end
   function _intM1i:getType() return "integrated" end
   function _intM1i:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1i = diags["M1i"].field
      specIn.volIntegral.compsVdim:advance(tm, {M1i}, {self.field})
   end

   -- ~~~~ Integrated particle energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2 = Proto(DiagsImplBase)
   function _intM2:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
      self:setVolIntegral(specIn)
   end
   function _intM2:getDependencies() return {"M2"} end
   function _intM2:getType() return "integrated" end
   function _intM2:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["M2"].field
      specIn.volIntegral.comps1:advance(tm, {M2}, {self.field})
   end

   -- ~~~~ Integrated mean flow energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Flow = Proto(DiagsImplBase)
   function _intM2Flow:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
      self:setVolIntegral(specIn)
   end
   function _intM2Flow:getDependencies() return {"M2Flow"} end
   function _intM2Flow:getType() return "integrated" end
   function _intM2Flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Flow = diags["M2Flow"].field
      specIn.volIntegral.comps1:advance(tm, {M2Flow}, {self.field})
   end

   -- ~~~~ Integrated thermal energy density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2Thermal = Proto(DiagsImplBase)
   function _intM2Thermal:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
      self:setVolIntegral(specIn)
   end
   function _intM2Thermal:getDependencies() return {"M2Thermal"} end
   function _intM2Thermal:getType() return "integrated" end
   function _intM2Thermal:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2Thermal = diags["M2Thermal"].field
      specIn.volIntegral.comps1:advance(tm, {M2Thermal}, {self.field})
   end

   -- ~~~~ L2 norm of the distribution function ~~~~~~~~~~~~~~~~~~~~~~
   local _intL2 = Proto(DiagsImplBase)
   function _intL2:fullInit(diagApp, specIn, field, owner)
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
      M0        = _M0,
      M1i       = _M1i,
      M2        = _M2,
      M2ij      = _M2ij,
      M3i       = _M3i,
      Udrift    = _Udrift,
      M2Flow    = _M2Flow,
      M2Thermal = _M2Thermal,
      VtSq      = _VtSq,
      intM0        = _intM0,
      intM1i       = _intM1i,
      intM2        = _intM2,
      intM2Flow    = _intM2Flow,
      intM2Thermal = _intM2Thermal,
      intL2        = _intL2,
   }
end

return implementation

