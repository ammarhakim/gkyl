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

-- The first entry is the diagnostic name. The second is the other diagnostics it depends on.
local allowedFieldDiags = {{"MomSq",{}}}
local allowedIntegratedDiags = {{"intMom",{}}}

local contains=function(tbl, elem) return lume.any(tbl, function(e) return e[1]==elem end) end

local FluidDiags = Proto(DiagsBase)  -- FluidDiags is a child of DiagsBase.

function FluidDiags:init()
   self.allowedFieldDiags = allowedFieldDiags
   self.allowedIntDiags   = allowedIntegratedDiags

   self.diagsImp = FluidDiags._diagImp
end
   
function FluidDiags:fullInit(mySpecies)

   self.name     = mySpecies.name
   self.grid     = mySpecies.grid
   self.basis    = mySpecies.basis
   self.nMoments = mySpecies.nMoments

   self.weakMultiply = mySpecies.weakMultiply
   self.weakDivide   = mySpecies.weakDivide

   self.fieldDiags, self.intDiags = {}, {}  -- Tables to store field and integrated diagnostics.

   -- Sort requested diagnostics into field and integrated diagnostics.
   for _, nm in ipairs(mySpecies.diagnostics) do
      if contains(self.allowedFieldDiags, nm) then
         self.fieldDiags[nm]      = self.diagsImp[nm]
         self.fieldDiags[nm].done = false
      elseif contains(self.allowedIntDiags, nm) then
         self.intDiags[nm]      = self.diagsImp[nm]
         self.intDiags[nm].done = false
      else
         assert(false, string.format("FluidDiagnostics: %s is not an allowed diagnostic.",nm))
      end
   end

   -- Organize diagnostics and identify dependencies.
   self:orgDiagnostics(self.allowedFieldDiags, self.fieldDiags)
   self:orgDiagnostics(self.allowedIntDiags, self.intDiags, self.allowedFieldDiags, self.fieldDiags)

   -- Initialize diagnostics.
   for _, diag in lume.orderedIter(self.fieldDiags) do diag.init(self, mySpecies) end
   for _, diag in lume.orderedIter(self.intDiags) do diag.init(self, mySpecies) end

   self.spentTime = { fieldDiags=0., intDiags=0. }
end

function FluidDiags:orgDiagnostics(allowed, diagsTbl, otherAllowed, otherDiagsTbl)
   -- Organize diagnostics. Insert other diagnostics needed. Set metatable for orderedIter.
   -- The last two arguments are intended to be the allowedField Diags and fieldDiags when
   -- the first two are the allowedIntegratedDiags and intDiags. That's because in such case
   -- organizing intDiags may change and reorder fieldDiags.
   assert((otherAllowed and otherDiagsTbl) or (otherAllowed==nil and otherDiagsTbl==nil),
          "GyrofluidDiags:orgDiagnostics: must pass both otherAllowed and otherDiagsTbl, or neither.")
   local reorderOther = false

   -- Check for dependencies and insert other necessary diagnostics.
   -- Note: this can result in extra diagnostics written out.
   for nm, diag in pairs(diagsTbl) do
      local info, _ = lume.match(allowed, function(e) return e[1]==nm end)
      local depends = info and info[2] or {}
      for depI, depNm in ipairs(depends) do
         if lume.find(diagsTbl, depNm)==nil then
            diagsTbl[depNm] = self.diagsImp[depNm]
         end
         -- Also check the other diagsTbl if present.
         if otherDiagsTbl then
            if lume.find(otherDiagsTbl, depNm)==nil then
               reorderOther = true
               otherDiagsTbl[depNm] = self.diagsImp[depNm]
            end
         end
      end
   end
   local function sortDiagsTbl(allowed, diags)
      -- Create a keys list so we always loop in the same order.
      local diags_keys = {}
      for nm in pairs(diags) do table.insert(diags_keys, nm) end
      local diagSort = function(d1,d2)
         -- Check if d1 is in the list of dependencies of d2. If so it must go before d2.
         local info, _  = lume.match(allowed, function(e) return e[1]==d2 end)
         local depends2 = info and info[2] or {}
         if lume.find(depends2, d1) then return true else return false end
      end
      table.sort(diags_keys, diagSort)
      setmetatable(diags, diags_keys)
   end
   sortDiagsTbl(allowed, diagsTbl)
   -- Also sort other diagsTbl if something was added.
   if otherDiagsTbl and reorderOther then sortDiagsTbl(otherAllowed, otherDiagsTbl) end
end

function FluidDiags:resetState(tm)
   -- Reset the state indicating if diagnostic has been computed.
   -- This state prevents computing a diagnostic twice when there are dependencies.
   for _, diag in lume.orderedIter(self.intDiags) do diag.done=false end
   for _, diag in lume.orderedIter(self.fieldDiags) do diag.done=false end
end

function FluidDiags:calcDiagDependencies(tm, mySpecies, diagNm)
   -- Given a diagnostic, compute the other diagnostics it may depend on.
   local info, _ = lume.match(self.allowedFieldDiags, function(e) return e[1]==diagNm end)
   local depends = info and info[2] or {}
   for _, depNm in ipairs(depends) do
      local depDiag = self.fieldDiags[depNm]
      if depDiag.done then break end
      self:calcDiagDependencies(tm, mySpecies, depNm)
      depDiag.calc(tm, mySpecies)
      depDiag.done = true
   end
end

function FluidDiags:calcIntegratedDiagnostics(tm, mySpecies)
   -- Compute integrated diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for diagNm, diag in lume.orderedIter(self.intDiags) do
      if diag.done then break end
      self:calcDiagDependencies(tm, spec, diagNm)
      diag.calc(tm, spec)
      diag.done = true
   end

   self.spentTime.intDiags = self.spentTime.intDiags + Time.clock() - tmStart
end

function FluidDiags:calcFieldDiagnostics(tm, mySpecies)
   -- Compute field diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for diagNm, diag in lume.orderedIter(self.fieldDiags) do
      if diag.done then break end
      self:calcDiagDependencies(tm, spec, diagNm)
      diag.calc(tm, spec)
      diag.done = true
   end

   self.spentTime.fieldDiags = self.spentTime.fieldDiags + Time.clock() - tmStart
end

function FluidDiags:writeFieldDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_%d.bp", self.name, nm, fr), tm, fr)
   end
end

function FluidDiags:writeIntegratedDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.intDiags) do
      diag.field:write(string.format("%s_%s.bp", self.name, nm), tm, fr)
   end
end

function FluidDiags:write(tm, fr)
   self:writeFieldDiagnostics(tm, fr)
   self:writeIntegratedDiagnostics(tm, fr)
end

function FluidDiags:writeRestartFieldDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_restart.bp", self.name, nm), tm, fr, false)
   end
end

function FluidDiags:writeRestartIntegratedDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.intDiags) do
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
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, nm))
   end

   for nm, diag in lume.orderedIter(self.intDiags) do
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
FluidDiags._diagImp["MomSq"].init = function(self, specIn)
   FluidDiags._diagImp["MomSq"].field    = specIn:allocMoment()
   FluidDiags._diagImp["MomSq"].updaters = specIn.weakMultiply
end
FluidDiags._diagImp["MomSq"].calc = function(tm, specIn)
   local momIn = specIn:getMoments(1)
   FluidDiags._diagImp["MomSq"].updaters:advance(tm, {momIn, momIn}, {FluidDiags._diagImp["MomSq"].field})
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
FluidDiags._diagImp["intMom"] = {}
FluidDiags._diagImp["intMom"].init = function(self, specIn)
   FluidDiags._diagImp["intMom"].field    = DataStruct.DynVector { numComponents = self.nMoments }
   FluidDiags._diagImp["intMom"].updaters = specIn.volIntegral.compsN
end
FluidDiags._diagImp["intMom"].calc = function(tm, specIn)
   FluidDiags._diagImp["intMom"].updaters:advance(tm, {specIn:rkStepperFields()[1]}, {FluidDiags._diagImp["intMom"].field})
end

return FluidDiags
