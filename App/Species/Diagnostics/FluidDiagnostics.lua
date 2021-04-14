-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for FluidSpecies.
-- 
-- Supported diagnostic are defined as functions at the bottom of the file.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local DiagsBase  = require "App.Species.Diagnostics.SpeciesDiagnosticsBase"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local lume       = require "Lib.lume"

local allowedFieldDiags = {"MomSq"}
local allowedIntegratedDiags = {"intMom"}

local function contains(tbl, elem) return lume.any(tbl, function(e) return e==elem end) end
local function isFieldDiagAllowed(nm) return contains(allowedFieldDiags, nm) end
local function isIntegratedDiagAllowed(nm) return contains(allowedIntegratedDiags, nm) end

local FluidDiags = Proto(DiagsBase)  -- FluidDiags is a child of DiagsBase.

function FluidDiags:init(mySpecies)

   self.name     = mySpecies.name
   self.grid     = mySpecies.grid
   self.basis    = mySpecies.basis
   self.nMoments = mySpecies.nMoments

   self.weakMultiply = mySpecies.weakMultiply
   self.weakDivide   = mySpecies.weakDivide
   
   self.fieldDiags, self.intDiags = {}, {}  -- Tables to store field and integrated diagnostics.

   -- Sort requested diagnostics into field and integrated diagnostics.
   for _, nm in ipairs(mySpecies.diagnostics) do
      if isFieldDiagAllowed(nm) then
         self.fieldDiags[nm] = FluidDiags._diagImp[nm]
         self.fieldDiags[nm].init(self, mySpecies)
      elseif isIntegratedDiagAllowed(nm) then
         self.intDiags[nm] = FluidDiags._diagImp[nm]
         self.intDiags[nm].init(self, mySpecies)
      else
         assert(false, string.format("FluidDiagnostics: %s is not an allowed diagnostic.",nm))
      end
   end

   self.spentTime = { fieldDiags=0., intDiags=0. }
end

function FluidDiags:calcIntegratedDiagnostics(tm, mySpecies)
   -- Compute integrated diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diag in pairs(self.intDiags) do diag.calc(tm, spec) end

   self.spentTime.intDiags = self.spentTime.intDiags + Time.clock() - tmStart
end

function FluidDiags:calcFieldDiagnostics(tm, mySpecies)
   -- Compute field diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diag in pairs(self.fieldDiags) do diag.calc(tm, spec) end

   self.spentTime.fieldDiags = self.spentTime.fieldDiags + Time.clock() - tmStart
end

function FluidDiags:writeFieldDiagnostics(tm, fr)
   for nm, diag in pairs(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_%d.bp", self.name, nm, fr), tm, fr)
   end
end

function FluidDiags:writeIntegratedDiagnostics(tm, fr)
   for nm, diag in pairs(self.intDiags) do
      diag.field:write(string.format("%s_%s.bp", self.name, nm), tm, fr)
   end
end

function FluidDiags:write(tm, fr)
   self:writeFieldDiagnostics(tm, fr)
   self:writeIntegratedDiagnostics(tm, fr)
end

function FluidDiags:writeRestartFieldDiagnostics(tm, fr)
   for nm, diag in pairs(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_restart.bp", self.name, nm), tm, fr, false)
   end
end

function FluidDiags:writeRestartIntegratedDiagnostics(tm, fr)
   for nm, diag in pairs(self.intDiags) do
      diag.field:write(string.format("%s_%s_restart.bp", self.name, nm), tm, fr, false, false)
   end
end

function FluidDiags:writeRestart(tm, fieldFr, intFr)
   self:writeRestartFieldDiagnostics(tm, fieldFr)

   -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
   -- restart write frequency is higher than the normal write frequency from nFrame.
   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self:writeRestartIntegratedDiagnostics(tm, intFr)
end

function FluidDiags:readRestart()
   local tm, fr
   for nm, diag in pairs(self.fieldDiags) do
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, nm))
   end

   for nm, diag in pairs(self.intDiags) do
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, nm))
   end
   return tm, fr
end

-----------------------------------
-- Place diagnostics below
--
FluidDiags._diagImp = {}

-- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
FluidDiags._diagImp["MomSq"] = {}
FluidDiags._diagImp["MomSq"].init = function(self, mySpecies)
   FluidDiags._diagImp["MomSq"].field    = mySpecies:allocMoment()
   FluidDiags._diagImp["MomSq"].updaters = Updater.CartFieldBinOp {
                                              onGrid    = self.grid,
                                              weakBasis = self.basis,
                                              operation = "Multiply",
                                           }
end
FluidDiags._diagImp["MomSq"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   FluidDiags._diagImp["MomSq"].updaters:advance(tm, {momIn, momIn}, {FluidDiags._diagImp["MomSq"].field})
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
FluidDiags._diagImp["intMom"] = {}
FluidDiags._diagImp["intMom"].init = function(self, mySpecies)
   FluidDiags._diagImp["intMom"].field    = DataStruct.DynVector { numComponents = self.nMoments }
   FluidDiags._diagImp["intMom"].updaters = Updater.CartFieldIntegratedQuantCalc {
                                               onGrid        = self.grid,
                                               basis         = self.basis,
                                               numComponents = self.nMoments,
                                               quantity      = "V"
                                            }
end
FluidDiags._diagImp["intMom"].calc = function(tm, specIn)
   FluidDiags._diagImp["intMom"].updaters:advance(tm, {specIn:rkStepperFields()[1]}, {FluidDiags._diagImp["intMom"].field})
end

return FluidDiags
