-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for GyrofluidSpecies.
-- 
-- We assume that integrated diagnostics can depend on field diagnostics, but
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

-- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_MomSq = Proto(DiagsImplBase)
function GyrofluidDiag_MomSq:fullInit(diagApp, specIn)
   self.field = specIn:allocVectorMoment(specIn.nMoments)
   self.done  = false
end
function GyrofluidDiag_MomSq:getType() return "grid" end
function GyrofluidDiag_MomSq:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getNoJacMoments()
   specIn.weakMultiply:advance(tm, {momIn, momIn}, {self.field})
end

-- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_M0 = Proto(DiagsImplBase)
function GyrofluidDiag_M0:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_M0:getType() return "grid" end
function GyrofluidDiag_M0:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getNoJacMoments()
   self.field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(1))
end

-- ~~~~ Momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_M1 = Proto(DiagsImplBase)
function GyrofluidDiag_M1:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_M1:getType() return "grid" end
function GyrofluidDiag_M1:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getNoJacMoments()
   self.field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(2))
end

-- ~~~~ Energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_M2 = Proto(DiagsImplBase)
function GyrofluidDiag_M2:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_M2:getType() return "grid" end
function GyrofluidDiag_M2:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getNoJacMoments()
   self.field:combineOffset(2./specIn:getMass(), momIn, specIn:getMomOff(3))
end

-- ~~~~ Flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_M2flow = Proto(DiagsImplBase)
function GyrofluidDiag_M2flow:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_M2flow:getDependencies() return {"M1","upar"} end
function GyrofluidDiag_M2flow:getType() return "grid" end
function GyrofluidDiag_M2flow:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local upar, M1      = diags["upar"].field, diags["M1"].field
   specIn.weakMultiply:advance(tm, {upar,M1}, {self.field})
end

-- ~~~~ Parallel flow speed ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_upar = Proto(DiagsImplBase)
function GyrofluidDiag_upar:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false

   self.fieldAux = {mJacM0=specIn:allocMoment(), mJacM1=specIn:allocMoment()}
end
function GyrofluidDiag_upar:getType() return "grid" end
function GyrofluidDiag_upar:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getMoments(1)
   specIn:uParCalc(tm, momIn, self.fieldAux.mJacM0, self.fieldAux.mJacM1, self.field)
end

-- ~~~~ Perpendicular pressure ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_pperp = Proto(DiagsImplBase)
function GyrofluidDiag_pperp:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false

   self.fieldAux = {M2perp=specIn:allocMoment()}
end
function GyrofluidDiag_pperp:getType() return "grid" end
function GyrofluidDiag_pperp:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getNoJacMoments()
   specIn:pPerpJacCalc(tm, momIn, self.fieldAux.M2perp, self.field)
end

-- ~~~~ Parallel pressure ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_ppar = Proto(DiagsImplBase)
function GyrofluidDiag_ppar:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false

   self.fieldAux = {mM2=specIn:allocMoment()}
end
function GyrofluidDiag_ppar:getDependencies() return {"pperp","M2flow"} end
function GyrofluidDiag_ppar:getType() return "grid" end
function GyrofluidDiag_ppar:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local pperp, M2flow = diags["pperp"].field, diags["M2flow"].field
   local momIn         = specIn:getNoJacMoments()
   self.fieldAux.mM2:combineOffset(1, momIn, specIn:getMomOff(3))
   self.field:combine( 2., self.fieldAux.mM2, -2., pperp, -specIn:getMass(), M2flow)
end

-- ~~~~ Perpendicular temperature ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_Tperp = Proto(DiagsImplBase)
function GyrofluidDiag_Tperp:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_Tperp:getDependencies() return {"M0","pperp"} end
function GyrofluidDiag_Tperp:getType() return "grid" end
function GyrofluidDiag_Tperp:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local pperp, M0     = diags["pperp"].field, diags["M0"].field
   specIn.weakDivide:advance(tm, {M0,pperp}, {self.field})
end

-- ~~~~ Parallel temperature ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_Tpar = Proto(DiagsImplBase)
function GyrofluidDiag_Tpar:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_Tpar:getDependencies() return {"M0","ppar"} end
function GyrofluidDiag_Tpar:getType() return "grid" end
function GyrofluidDiag_Tpar:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local ppar, M0      = diags["ppar"].field, diags["M0"].field
   specIn.weakDivide:advance(tm, {M0,ppar}, {self.field})
end

-- ~~~~ Scalar temperature ~~~~~~~~~~~~~~~~~~~~~~
-- MF 2021/04/14: I think it is more efficient to implement this using
--                3*n*T = 2*E - m*n*upar^2.
local GyrofluidDiag_Temp = Proto(DiagsImplBase)
function GyrofluidDiag_Temp:fullInit(diagApp, specIn)
   self.field = specIn:allocMoment()
   self.done  = false
end
function GyrofluidDiag_Temp:getDependencies() return {"M0","M2flow"} end
function GyrofluidDiag_Temp:getType() return "grid" end
function GyrofluidDiag_Temp:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local M2flow, M0    = diags["M2flow"].field, diags["M0"].field
   self.field:combineOffset( 2./3., specIn:getNoJacMoments(), specIn:getMomOff(3),
                            -specIn:getMass()/3., M2flow, 0)
   specIn.weakDivide:advance(tm, {M0,self.field}, {self.field})
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_intMom = Proto(DiagsImplBase)
function GyrofluidDiag_intMom:fullInit(diagApp, specIn)
   self.field    = DataStruct.DynVector { numComponents = specIn.nMoments }
   self.updaters = specIn.volIntegral.compsN
   self.done     = false
end
function GyrofluidDiag_intMom:getType() return "integrated" end
function GyrofluidDiag_intMom:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   self.updaters:advance(tm, {specIn:getNoJacMoments()}, {self.field})
end

-- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_intM0 = Proto(DiagsImplBase)
function GyrofluidDiag_intM0:fullInit(diagApp, specIn)
   self.field = DataStruct.DynVector { numComponents = 1 }
   self.done  = false
end
function GyrofluidDiag_intM0:getDependencies() return {"M0"} end
function GyrofluidDiag_intM0:getType() return "integrated" end
function GyrofluidDiag_intM0:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local M0 = diags["M0"].field
   specIn.volIntegral.comps1:advance(tm, {M0}, {self.field})
end

-- ~~~~ Integrated momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_intM1 = Proto(DiagsImplBase)
function GyrofluidDiag_intM1:fullInit(diagApp, specIn)
   self.field = DataStruct.DynVector { numComponents = 1 }
   self.done  = false
end
function GyrofluidDiag_intM1:getDependencies() return {"M1"} end
function GyrofluidDiag_intM1:getType() return "integrated" end
function GyrofluidDiag_intM1:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local M1 = diags["M1"].field
   specIn.volIntegral.comps1:advance(tm, {M1}, {self.field})
end

-- ~~~~ Integrated energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_intM2 = Proto(DiagsImplBase)
function GyrofluidDiag_intM2:fullInit(diagApp, specIn)
   self.field = DataStruct.DynVector { numComponents = 1 }
   self.done  = false
end
function GyrofluidDiag_intM2:getDependencies() return {"M2"} end
function GyrofluidDiag_intM2:getType() return "integrated" end
function GyrofluidDiag_intM2:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local M2 = diags["M2"].field
   specIn.volIntegral.comps1:advance(tm, {M2}, {self.field})
end

-- ~~~~ Integrated flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
local GyrofluidDiag_intM2flow = Proto(DiagsImplBase)
function GyrofluidDiag_intM2flow:fullInit(diagApp, specIn)
   self.field = DataStruct.DynVector { numComponents = 1 }
   self.done  = false
end
function GyrofluidDiag_intM2flow:getDependencies() return {"M2flow"} end
function GyrofluidDiag_intM2flow:getType() return "integrated" end
function GyrofluidDiag_intM2flow:advance(tm, inFlds, outFlds)
   local specIn, diags = inFlds[1], inFlds[2]
   local M2flow = diags["M2flow"].field
   specIn.volIntegral.comps1:advance(tm, {M2flow}, {self.field})
end

return {
   MomSq  = GyrofluidDiag_MomSq,
   M0     = GyrofluidDiag_M0,
   M1     = GyrofluidDiag_M1,
   M2     = GyrofluidDiag_M2,
   M2flow = GyrofluidDiag_M2flow,
   upar   = GyrofluidDiag_upar,
   pperp  = GyrofluidDiag_pperp,
   ppar   = GyrofluidDiag_ppar,
   Tperp  = GyrofluidDiag_Tperp,
   Tpar   = GyrofluidDiag_Tpar,
   Temp   = GyrofluidDiag_Temp,
   intMom    = GyrofluidDiag_intMom,
   intM0     = GyrofluidDiag_intM0,
   intM1     = GyrofluidDiag_intM1,
   intM2     = GyrofluidDiag_intM2,
   intM2flow = GyrofluidDiag_intM2flow,
}

