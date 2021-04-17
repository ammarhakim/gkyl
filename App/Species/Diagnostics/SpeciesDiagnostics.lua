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
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local lume       = require "Lib.lume"

local contains=function(tbl, elem)
   for k, _ in pairs(tbl) do if elem == k then return true end end
   return false
end

local SpeciesDiagnostics = Proto()  -- FluidDiags is a child of DiagsBase.

function SpeciesDiagnostics:init(tbl)
   self.tbl = tbl
end

function SpeciesDiagnostics:fullInit(mySpecies, diagsImpl)

   self.name = mySpecies.name

   self.diags = {}  -- Grid and integrated diagnostics.
   self.diagGroups = {grid={}, integrated={}}  -- Names of requested diags of each kind.

   -- Sort requested diagnostics into grid and integrated diagnostics.
   for _, nm in ipairs(mySpecies.tbl.diagnostics) do
      if contains(diagsImpl, nm) then
         self.diags[nm] = diagsImpl[nm]{}
         local diagType = diagsImpl[nm]:getType()
         table.insert(self.diagGroups[diagType], nm)
      else
         assert(false, string.format("FluidDiagnostics: %s is not an allowed diagnostic.",nm))
      end
   end

   -- Organize diagnostics and identify dependencies.
   self:orgDiagnostics(diagsImpl, self.diags)

   -- Initialize diagnostics.
   for _, dG in pairs(self.diagGroups) do
      for _, diagNm in ipairs(dG) do self.diags[diagNm]:fullInit(self, mySpecies) end
   end

   self.spentTime = { gridDiags=0., intDiags=0. }
end

function SpeciesDiagnostics:orgDiagnostics(diagsImpl, diagsTbl)
   -- Insert other diagnostics needed by the requested diagnostics. Then organize them
   -- so the dependents are computed last. 

   -- Check for dependencies and insert other necessary diagnostics.
   -- Note: this can result in extra diagnostics written out.
   local function checkNinsertDependencies(diagNm, impl, diags)
      local depends = impl[diagNm]:getDependencies()
      -- Loop over dependencies and check if they are in table. If not, include them.
      for depI, depNm in ipairs(depends) do
         if diags[depNm]==nil then
            if contains(impl, depNm) then
               diags[depNm] = impl[depNm]
               table.insert(self.diagGroups[impl[depNm]:getType()], depNm) 
               checkNinsertDependencies(depNm, impl, diags)
            end
         end
      end
   end
   for nm, _ in pairs(diagsTbl) do checkNinsertDependencies(nm, diagsImpl, diagsTbl) end

   local function sortDiagsTbl(impl, diagList)
      -- Create a keys list so we always loop in the same order.
      local sortFunc = function(d1,d2)
         -- Check if d1 is in the list of dependencies of d2, or the dependencies of
         -- the dependencies of d2, and so on (recursion). If so it must go before d2.
         local function dependsOn(diag1,diag2)
            local depends = impl[diag2]:getDependencies()
            local found   = false  
            if lume.find(depends, diag1) then found = true end
            if not found then
               for _, oNm in ipairs(depends) do found = found or dependsOn(diag1, oNm) end
            end
            return found
         end
         return dependsOn(d1,d2)
      end
      table.sort(diagList, sortFunc)
   end
   sortDiagsTbl(diagsImpl, self.diagGroups.grid)
   sortDiagsTbl(diagsImpl, self.diagGroups.integrated)
end

function SpeciesDiagnostics:resetState(tm)
   -- Reset the state indicating if diagnostic has been computed.
   -- This state prevents computing a diagnostic twice when there are dependencies.
   for _, diagNm in ipairs(self.diagGroups.grid) do self.diags[diagNm].done=false end
   for _, diagNm in ipairs(self.diagGroups.integrated) do self.diags[diagNm].done=false end
end

function SpeciesDiagnostics:calcDiagDependencies(tm, mySpecies, diagNm)
   -- Given a diagnostic, compute the other diagnostics it may depend on.
   local spec    = mySpecies
   local depends = self.diags[diagNm]:getDependencies()

   for _, depNm in ipairs(depends) do
      local depDiag = self.diags[depNm]
      if not depDiag.done then
         self:calcDiagDependencies(tm, spec, depNm)
         depDiag:advance(tm, {spec, self.diags}, {})
         depDiag.done = true
      end
   end
end

function SpeciesDiagnostics:calcIntegratedDiagnostics(tm, mySpecies)
   -- Compute integrated diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      if not diag.done then
         self:calcDiagDependencies(tm, spec, diagNm)
         diag:advance(tm, {spec, self.diags}, {})
         diag.done = true
      end
   end

   self.spentTime.intDiags = self.spentTime.intDiags + Time.clock() - tmStart
end

function SpeciesDiagnostics:calcGridDiagnostics(tm, mySpecies)
   -- Compute grid diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      if not diag.done then
         self:calcDiagDependencies(tm, spec, diagNm)
         diag:advance(tm, {spec, self.diags}, {})
         diag.done = true
      end
   end

   self.spentTime.gridDiags = self.spentTime.gridDiags + Time.clock() - tmStart
end

function SpeciesDiagnostics:writeGridDiagnostics(tm, fr)
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s_%d.bp", self.name, diagNm, fr), tm, fr)
   end
end

function SpeciesDiagnostics:writeIntegratedDiagnostics(tm, fr)
   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s.bp", self.name, diagNm), tm, fr)
   end
end

function SpeciesDiagnostics:write(tm, fr)
   self:writeGridDiagnostics(tm, fr)
   self:writeIntegratedDiagnostics(tm, fr)
end

function SpeciesDiagnostics:writeRestartGridDiagnostics(tm, fr)
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s_restart.bp", self.name, diagNm), tm, fr, false)
   end
end

function SpeciesDiagnostics:writeRestartIntegratedDiagnostics(tm, fr)
   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s_restart.bp", self.name, diagNm), tm, fr, false, false)
   end
end

function SpeciesDiagnostics:writeRestart(tm, gridFr, intFr)
   self:writeRestartGridDiagnostics(tm, gridFr)

   -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
   -- restart write frequency is higher than the normal write frequency from nFrame.
   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self:writeRestartIntegratedDiagnostics(tm, intFr)
end

function SpeciesDiagnostics:readRestart()
   local tm, fr
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, diagNm))
   end

   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, diagNm))
   end
   return tm, fr
end

function SpeciesDiagnostics:getDiagTime()
   return self.spentTime.gridDiags + self.spentTime.intDiags
end

return SpeciesDiagnostics
