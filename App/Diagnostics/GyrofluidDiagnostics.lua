-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for GyrofluidSpecies.
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
   -- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
   local _MomSq = Proto(DiagsImplBase)
   function _MomSq:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocVectorMoment(specIn.nMoments)
      self.owner = owner
      self.done  = false
   end
   function _MomSq:getType() return "grid" end
   function _MomSq:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = self.owner:getNoJacMoments()
      specIn.weakMultiply:advance(tm, {momIn, momIn}, {self.field})
   end
   
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
      local momIn  = self.owner:getNoJacMoments()
      self.field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(1))
   end
   
   -- ~~~~ Momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
   local _M1 = Proto(DiagsImplBase)
   function _M1:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M1:getType() return "grid" end
   function _M1:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = self.owner:getNoJacMoments()
      self.field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(2))
   end
   
   -- ~~~~ Energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
   local _M2 = Proto(DiagsImplBase)
   function _M2:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M2:getType() return "grid" end
   function _M2:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = self.owner:getNoJacMoments()
      self.field:combineOffset(2./specIn:getMass(), momIn, specIn:getMomOff(3))
   end
   
   -- ~~~~ Flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
   local _M2flow = Proto(DiagsImplBase)
   function _M2flow:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _M2flow:getDependencies() return {"M1","upar"} end
   function _M2flow:getType() return "grid" end
   function _M2flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local upar, M1      = diags["upar"].field, diags["M1"].field
      specIn.weakMultiply:advance(tm, {upar,M1}, {self.field})
   end
   
   -- ~~~~ Parallel flow speed ~~~~~~~~~~~~~~~~~~~~~~
   local _upar = Proto(DiagsImplBase)
   function _upar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   
      self.fieldAux = {mJacM0=specIn:allocMoment(), mJacM1=specIn:allocMoment()}
   end
   function _upar:getType() return "grid" end
   function _upar:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = specIn:getMoments(1)
      specIn:uParCalc(tm, momIn, self.fieldAux.mJacM0, self.fieldAux.mJacM1, self.field)
   end
   
   -- ~~~~ Perpendicular pressure ~~~~~~~~~~~~~~~~~~~~~~
   local _pperp = Proto(DiagsImplBase)
   function _pperp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   
      self.fieldAux = {M2perp=specIn:allocMoment()}
   end
   function _pperp:getType() return "grid" end
   function _pperp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = self.owner:getNoJacMoments()
      specIn:pPerpJacCalc(tm, momIn, self.fieldAux.M2perp, self.field)
   end
   
   -- ~~~~ Parallel pressure ~~~~~~~~~~~~~~~~~~~~~~
   local _ppar = Proto(DiagsImplBase)
   function _ppar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   
      self.fieldAux = {mM2=specIn:allocMoment()}
   end
   function _ppar:getDependencies() return {"pperp","M2flow"} end
   function _ppar:getType() return "grid" end
   function _ppar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local pperp, M2flow = diags["pperp"].field, diags["M2flow"].field
      local momIn         = self.owner:getNoJacMoments()
      self.fieldAux.mM2:combineOffset(1, momIn, specIn:getMomOff(3))
      self.field:combine( 2., self.fieldAux.mM2, -2., pperp, -specIn:getMass(), M2flow)
   end
   
   -- ~~~~ Perpendicular temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tperp = Proto(DiagsImplBase)
   function _Tperp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _Tperp:getDependencies() return {"M0","pperp"} end
   function _Tperp:getType() return "grid" end
   function _Tperp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local pperp, M0     = diags["pperp"].field, diags["M0"].field
      specIn.weakDivide:advance(tm, {M0,pperp}, {self.field})
   end
   
   -- ~~~~ Parallel temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tpar = Proto(DiagsImplBase)
   function _Tpar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _Tpar:getDependencies() return {"M0","ppar"} end
   function _Tpar:getType() return "grid" end
   function _Tpar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local ppar, M0      = diags["ppar"].field, diags["M0"].field
      specIn.weakDivide:advance(tm, {M0,ppar}, {self.field})
   end
   
   -- ~~~~ Scalar temperature ~~~~~~~~~~~~~~~~~~~~~~
   -- MF 2021/04/14: If Tpar and Tperp are not requested, I think it is more
   --                efficient to implement this using 3*n*T = 2*E - m*n*upar^2.
   local _Temp = Proto(DiagsImplBase)
   function _Temp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _Temp:getDependencies() return {"M0","M2flow"} end
   function _Temp:getType() return "grid" end
   function _Temp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2flow, M0    = diags["M2flow"].field, diags["M0"].field
      self.field:combineOffset( 2./3., self.owner:getNoJacMoments(), specIn:getMomOff(3),
                               -specIn:getMass()/3., M2flow, 0)
      specIn.weakDivide:advance(tm, {M0,self.field}, {self.field})
   end
   
   -- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
   local _intMom = Proto(DiagsImplBase)
   function _intMom:fullInit(diagApp, specIn, field, owner)
      self.field    = DataStruct.DynVector { numComponents = specIn.nMoments }
      self.updaters = specIn.volIntegral.compsN
      self.owner    = owner
      self.done     = false
   end
   function _intMom:getType() return "integrated" end
   function _intMom:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      self.updaters:advance(tm, {self.owner:getNoJacMoments()}, {self.field})
   end
   
   -- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM0:getDependencies() return {"M0"} end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0 = diags["M0"].field
      specIn.volIntegral.comps1:advance(tm, {M0}, {self.field})
   end
   
   -- ~~~~ Integrated momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
   local _intM1 = Proto(DiagsImplBase)
   function _intM1:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM1:getDependencies() return {"M1"} end
   function _intM1:getType() return "integrated" end
   function _intM1:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M1 = diags["M1"].field
      specIn.volIntegral.comps1:advance(tm, {M1}, {self.field})
   end
   
   -- ~~~~ Integrated energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2 = Proto(DiagsImplBase)
   function _intM2:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM2:getDependencies() return {"M2"} end
   function _intM2:getType() return "integrated" end
   function _intM2:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2 = diags["M2"].field
      specIn.volIntegral.comps1:advance(tm, {M2}, {self.field})
   end
   
   -- ~~~~ Integrated flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
   local _intM2flow = Proto(DiagsImplBase)
   function _intM2flow:fullInit(diagApp, specIn, field, owner)
      self.field = DataStruct.DynVector { numComponents = 1 }
      self.done  = false
   end
   function _intM2flow:getDependencies() return {"M2flow"} end
   function _intM2flow:getType() return "integrated" end
   function _intM2flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2flow = diags["M2flow"].field
      specIn.volIntegral.comps1:advance(tm, {M2flow}, {self.field})
   end
   
   return {
      MomSq  = _MomSq,
      M0     = _M0,
      M1     = _M1,
      M2     = _M2,
      M2flow = _M2flow,
      upar   = _upar,
      pperp  = _pperp,
      ppar   = _ppar,
      Tperp  = _Tperp,
      Tpar   = _Tpar,
      Temp   = _Temp,
      intMom    = _intMom,
      intM0     = _intM0,
      intM1     = _intM1,
      intM2     = _intM2,
      intM2flow = _intM2flow,
   }
end

return implementation
