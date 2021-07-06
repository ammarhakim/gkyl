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
   function _M2flow:getDependencies() return {"M1","Upar"} end
   function _M2flow:getType() return "grid" end
   function _M2flow:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local Upar, M1      = diags["Upar"].field, diags["M1"].field
      specIn.weakMultiply:advance(tm, {Upar,M1}, {self.field})
   end
   
   -- ~~~~ Perpendicular energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
   local _M2perp = Proto(DiagsImplBase)
   function _M2perp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _M2perp:getDependencies() return {"Pperp"} end
   function _M2perp:getType() return "grid" end
   function _M2perp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local Pperp = diags["Pperp"].field
      self.field:combine(2./specIn:getMass(), Pperp)
   end
   
   -- ~~~~ Parallel energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
   local _M2par = Proto(DiagsImplBase)
   function _M2par:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _M2par:getDependencies() return {"M2","M2perp"} end
   function _M2par:getType() return "grid" end
   function _M2par:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M2, M2perp = diags["M2"].field, diags["M2perp"].field
      self.field:combine(1., M2, -1., M2perp)
   end
   
   -- ~~~~ Parallel flow speed ~~~~~~~~~~~~~~~~~~~~~~
   local _Upar = Proto(DiagsImplBase)
   function _Upar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   
      self.fieldAux = {mJacM0=specIn:allocMoment(), mJacM1=specIn:allocMoment()}
   end
   function _Upar:getType() return "grid" end
   function _Upar:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = specIn:getMoments(1)
      specIn:uParCalc(tm, momIn, self.fieldAux.mJacM0, self.fieldAux.mJacM1, self.field)
   end
   
   -- ~~~~ Perpendicular pressure ~~~~~~~~~~~~~~~~~~~~~~
   local _Pperp = Proto(DiagsImplBase)
   function _Pperp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   
      self.fieldAux = {M2perp=specIn:allocMoment()}
   end
   function _Pperp:getType() return "grid" end
   function _Pperp:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = self.owner:getNoJacMoments()
      specIn:pPerpJacCalc(tm, momIn, self.fieldAux.M2perp, self.field)
   end
   
   -- ~~~~ Parallel pressure ~~~~~~~~~~~~~~~~~~~~~~
   local _Ppar = Proto(DiagsImplBase)
   function _Ppar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.owner = owner
      self.done  = false
   
      self.fieldAux = {mM2=specIn:allocMoment()}
   end
   function _Ppar:getDependencies() return {"Pperp","M2flow"} end
   function _Ppar:getType() return "grid" end
   function _Ppar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local Pperp, M2flow = diags["Pperp"].field, diags["M2flow"].field
      local momIn         = self.owner:getNoJacMoments()
      self.fieldAux.mM2:combineOffset(1, momIn, specIn:getMomOff(3))
      self.field:combine( 2., self.fieldAux.mM2, -2., Pperp, -specIn:getMass(), M2flow)
   end
   
   -- ~~~~ Perpendicular temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tperp = Proto(DiagsImplBase)
   function _Tperp:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _Tperp:getDependencies() return {"M0","Pperp"} end
   function _Tperp:getType() return "grid" end
   function _Tperp:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local Pperp, M0     = diags["Pperp"].field, diags["M0"].field
      specIn.weakDivide:advance(tm, {M0,Pperp}, {self.field})
   end
   
   -- ~~~~ Parallel temperature ~~~~~~~~~~~~~~~~~~~~~~
   local _Tpar = Proto(DiagsImplBase)
   function _Tpar:fullInit(diagApp, specIn, field, owner)
      self.field = specIn:allocMoment()
      self.done  = false
   end
   function _Tpar:getDependencies() return {"M0","Ppar"} end
   function _Tpar:getType() return "grid" end
   function _Tpar:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local Ppar, M0      = diags["Ppar"].field, diags["M0"].field
      specIn.weakDivide:advance(tm, {M0,Ppar}, {self.field})
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
      self.updaters = specIn.volIntegral.vector
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
      specIn.volIntegral.scalar:advance(tm, {M0}, {self.field})
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
      specIn.volIntegral.scalar:advance(tm, {M1}, {self.field})
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
      specIn.volIntegral.scalar:advance(tm, {M2}, {self.field})
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
      specIn.volIntegral.scalar:advance(tm, {M2flow}, {self.field})
   end
   
   return {
      MomSq  = _MomSq,
      M0     = _M0,
      M1     = _M1,
      M2     = _M2,
      M2flow = _M2flow,
      M2perp = _M2perp,
      M2par  = _M2par,
      Upar   = _Upar,
      Pperp  = _Pperp,
      Ppar   = _Ppar,
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
