-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for FluidSpecies.
-- 
-- We assume that integrated diagnostics can depend on field diagnostics, but
-- not the other way around. We also assume that the integrated diagnostics are
-- always computed before the field diagnostics, since they may be computed more
-- frequently. This allows us to reset the state in calcIntegratedDiagostics.
-- 
-- Supported diagnostic are defined as functions at the bottom of the file.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local FluidDiags = require "App.Species.Diagnostics.FluidDiagnostics"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local lume       = require "Lib.lume"

-- The first entry is the diagnostic name. The second is the other diagnostics it depends on.
local allowedFieldDiags = {{"MomSq",{}},{"M0",{}},{"M1",{}},{"M2",{}}, {"M2flow",{"M1","upar"}},
                           {"upar",{}}, {"Tpar",{}}, {"Tperp",{}}, {"Temp",{}}, {"ppar",{}}, {"pperp",{}}}
local allowedIntegratedDiags = {{"intMom",{}},{"intM0",{"M0"}}, {"intM1",{"M1"}}, {"intM2",{"M2"}}, {"intM2flow",{"M2flow"}}}

local GyrofluidDiags = Proto(FluidDiags)  -- GyrofluidDiags is a child of FluidDiagnostics.

function GyrofluidDiags:init()
   self.allowedFieldDiags = allowedFieldDiags
   self.allowedIntDiags   = allowedIntegratedDiags

   self.diagsImp = GyrofluidDiags._diagImp
end
   
-----------------------------------
-- Place diagnostics below
--
GyrofluidDiags._diagImp = {}

-- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["MomSq"] = {}
GyrofluidDiags._diagImp["MomSq"].init = function(self, specIn)
   GyrofluidDiags._diagImp["MomSq"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["MomSq"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   specIn.weakMultiply:advance(tm, {momIn, momIn}, {GyrofluidDiags._diagImp["MomSq"].field})
end

-- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M0"] = {}
GyrofluidDiags._diagImp["M0"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M0"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M0"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M0"].field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(1))
end

-- ~~~~ Momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M1"] = {}
GyrofluidDiags._diagImp["M1"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M1"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M1"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M1"].field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(2))
end

-- ~~~~ Energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M2"] = {}
GyrofluidDiags._diagImp["M2"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M2"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M2"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M2"].field:combineOffset(2./specIn:getMass(), momIn, specIn:getMomOff(3))
end

-- ~~~~ Flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M2flow"] = {}
GyrofluidDiags._diagImp["M2flow"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M2flow"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M2flow"].calc = function(tm, specIn)
   specIn.weakMultiply:advance(tm, {GyrofluidDiags._diagImp["upar"].field,GyrofluidDiags._diagImp["M1"].field},
                               {GyrofluidDiags._diagImp["M2flow"].field})
end

-- ~~~~ Parallel flow speed ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["upar"] = {}
GyrofluidDiags._diagImp["upar"].init = function(self, specIn)
   GyrofluidDiags._diagImp["upar"].field = specIn:allocMoment()
   GyrofluidDiags._diagImp["upar"].fieldAux = {mJacM0=specIn:allocMoment(),mJacM1=specIn:allocMoment()}
end
GyrofluidDiags._diagImp["upar"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   specIn:uParCalc(tm, momIn, GyrofluidDiags._diagImp["upar"].fieldAux.mJacM0,
                              GyrofluidDiags._diagImp["upar"].fieldAux.mJacM1,
                              GyrofluidDiags._diagImp["upar"].field)
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["intMom"] = {}
GyrofluidDiags._diagImp["intMom"].init = function(self, specIn)
   GyrofluidDiags._diagImp["intMom"].field    = DataStruct.DynVector { numComponents = self.nMoments }
   GyrofluidDiags._diagImp["intMom"].updaters = specIn.volIntegral.compsN
end
GyrofluidDiags._diagImp["intMom"].calc = function(tm, specIn)
   GyrofluidDiags._diagImp["intMom"].updaters:advance(tm, {specIn:rkStepperFields()[1]}, {GyrofluidDiags._diagImp["intMom"].field})
end

-- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["intM0"] = {}
GyrofluidDiags._diagImp["intM0"].init = function(self, specIn)
   GyrofluidDiags._diagImp["intM0"].field    = DataStruct.DynVector { numComponents = 1 }
   GyrofluidDiags._diagImp["intM0"].updaters = specIn.volIntegral.comps1
end
GyrofluidDiags._diagImp["intM0"].calc = function(tm, specIn)
   GyrofluidDiags._diagImp["intM0"].updaters:advance(tm, {GyrofluidDiags._diagImp["M0"].field}, {GyrofluidDiags._diagImp["intM0"].field})
end

return GyrofluidDiags
