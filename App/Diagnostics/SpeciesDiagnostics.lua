-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for any species, or other objects in a species.
-- 
-- Supported diagnostics are defined elsewhere.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local lume       = require "Lib.lume"
local Mpi        = require "Comm.Mpi"
local ffi        = require "ffi"
local xsys       = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

local new, copy = xsys.from(ffi,"new, copy")

local contains=function(tbl, elem)
   for k, _ in pairs(tbl) do if elem == k then return true end end
   return false
end

local function orgDiagnostics(comm, diagsImpl, diagsTbl, diagGroups)
   -- Insert other diagnostics needed by the requested diagnostics. Then organize them
   -- so the dependents are computed last. 

   -- Check for dependencies and insert other necessary diagnostics.
   -- Note: this can result in extra diagnostics written out.
   local function checkNinsertDependencies(diagNm, impl, diags)
      local depends = impl[diagNm]:getDependencies()
      -- Loop over dependencies and check if they are in table. If not, include them.
      for _, depNm in ipairs(depends) do
         if diags[depNm]==nil then
            if contains(impl, depNm) then
               diags[depNm] = impl[depNm]
               table.insert(diagGroups[impl[depNm]:getType()], depNm) 
               checkNinsertDependencies(depNm, impl, diags)
            end
         end
      end
   end

   local diagsCurr = {}  -- Need a proxy table to create the iterator below, since diagsTbl will change.
   for nm, _ in pairs(diagsTbl) do table.insert(diagsCurr, nm) end
   for _, nm in ipairs(diagsCurr) do checkNinsertDependencies(nm, diagsImpl, diagsTbl) end

   -- Create a table of unique IDs. This is necessary to ensure sorting (below)
   -- has the same result in all MPI ranks. We tried broadcasting a string of the
   -- sorted concatenated diag names but that was not functioning reliably.
   local diag_keys = {}
   for k in pairs(diagsImpl) do table.insert(diag_keys, k) end
   table.sort(diag_keys)
   local diagIDs = {}
   for id, nm in ipairs(diag_keys) do diagIDs[nm] = id end

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
      -- Different MPI ranks may find different acceptable orders.
      -- Need to sync this order across MPI ranks.
      local myOrderedIDs = {}
      for _, nm in ipairs(diagList) do table.insert(myOrderedIDs, diagIDs[nm]) end
      local numIDs = #myOrderedIDs
      if numIDs > 0 then
         local Cint = new("int [?]", numIDs)
         for i, id in ipairs(myOrderedIDs) do Cint[i-1] = id end
--         Mpi.Bcast(Cint, numIDs, Mpi.INT, 0, Mpi.COMM_WORLD)
         Mpi.Bcast(Cint, numIDs, Mpi.INT, 0, comm)
         for i = 1,numIDs do myOrderedIDs[i] = Cint[i-1] end
         for i, id in ipairs(myOrderedIDs) do diagList[i] = lume.find(diagIDs, id) end 
      end
   end
   -- We have to sort the groups in the same order in every MPI rank.
   for _, v in lume.orderedIter(diagGroups) do sortDiagsTbl(diagsImpl, v) end
end

local SpeciesDiagnostics = Proto()

function SpeciesDiagnostics:init(tbl) self.tbl = tbl end

function SpeciesDiagnostics:fullInit(mySpecies, field, diagOwner)

   local diagsImpl = assert(self.tbl.implementation,
      "SpeciesDiagnostics: must specify the implementation of the diagnostics in 'implementation'.")

   self.name  = diagOwner.name

   self.diags = {}  -- Grid and integrated diagnostics.
   self.diagGroups = {grid={}, integrated={}}  -- Names of requested diags of each kind.
   lume.setOrder(self.diagGroups)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   -- Option to write grid and integrated diagnostics in two files, instead of every diagnostic having a file.
   self.inTwoFiles = false
   if lume.any(diagOwner.tbl.diagnostics, function(e) return e=="twoFiles" end) or
      lume.any(diagOwner.tbl.diagnostics, function(e) return e=="groupDiagnostics" end) or
      mySpecies.groupDiags then
      lume.remove(diagOwner.tbl.diagnostics,"twoFiles")
      lume.remove(diagOwner.tbl.diagnostics,"groupDiagnostics")
      self.inTwoFiles = true
   end

   -- Sort requested diagnostics into grid and integrated diagnostics.   
   for _, nm in ipairs(diagOwner.tbl.diagnostics) do
      if contains(diagsImpl, nm) then
         self.diags[nm] = diagsImpl[nm]{}
         local diagType = diagsImpl[nm]:getType()
         table.insert(self.diagGroups[diagType], nm)
      else
         assert(false, string.format("SpeciesDiagnostics: %s is not an allowed diagnostic.",nm))
      end
   end

   -- Organize diagnostics and identify dependencies.
   local comm = mySpecies.confGrid:commSet().comm
   orgDiagnostics(comm, diagsImpl, self.diags, self.diagGroups)

   -- Initialize diagnostics.
   for _, dG in pairs(self.diagGroups) do
      for _, diagNm in ipairs(dG) do self.diags[diagNm]:fullInit(self, mySpecies, field, diagOwner) end
   end

   self.writeGhost = xsys.pickBool(diagOwner.writeGhost, false)

   -- Pre-define methods for writing diagnostics to two files only (grid and integrated)
   -- or to write each diagnostic to its own file.
   -- MF 2021/06/06: Writing all integrated diagnostics to a single file is currently not available.
   if self.inTwoFiles then
      self.gridDiagsToWrite, self.gridDiagsToRead = {}, {}
      self.ioMethod = "MPI"
      local elemType
      for _, diagNm in ipairs(self.diagGroups["grid"]) do
         local diag = self.diags[diagNm]
         local fldElemType = diag.field:elemType()
         elemType = elemType or fldElemType
         assert(elemType==fldElemType, "SpeciesDiagnostics: elemType of all diagnostics must be the same.")
      end
      self.diagIo = AdiosCartFieldIo {
         elemType   = elemType,
         method     = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData   = {polyOrder = mySpecies.basis:polyOrder(),
                       basisType = mySpecies.basis:id(),},
      }
      self.writeFunc = function(tm, fr)
         SpeciesDiagnostics["writeGridDiagnostics_oneFile"](self, tm, fr)
         SpeciesDiagnostics["writeIntegratedDiagnostics_separateFiles"](self, tm, fr)
      end
      self.writeRestartFunc = function(tm, gridFr, intFr)
         SpeciesDiagnostics["writeRestartGridDiagnostics_oneFile"](self, tm, gridFr)

         -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
         -- restart write frequency is higher than the normal write frequency from nFrame.
         -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
         SpeciesDiagnostics["writeRestartIntegratedDiagnostics_separateFiles"](self, tm, intFr)
      end
      self.readRestartFunc = function()
         local tm, fr = SpeciesDiagnostics["readRestartGridDiagnostics_oneFile"](self)
         SpeciesDiagnostics["readRestartIntegratedDiagnostics_separateFiles"](self)
         return tm, fr
      end
   else
      self.writeFunc = function(tm, fr)
         SpeciesDiagnostics["writeGridDiagnostics_separateFiles"](self, tm, fr)
         SpeciesDiagnostics["writeIntegratedDiagnostics_separateFiles"](self, tm, fr)
      end
      self.writeRestartFunc = function(tm, gridFr, intFr)
         SpeciesDiagnostics["writeRestartGridDiagnostics_separateFiles"](self, tm, gridFr)

         -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
         -- restart write frequency is higher than the normal write frequency from nFrame.
         -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
         SpeciesDiagnostics["writeRestartIntegratedDiagnostics_separateFiles"](self, tm, intFr)
      end
      self.readRestartFunc = function()
         local tm, fr = SpeciesDiagnostics["readRestartGridDiagnostics_separateFiles"](self)
         SpeciesDiagnostics["readRestartIntegratedDiagnostics_separateFiles"](self)
         return tm, fr
      end
   end

   self.spentTime = { gridDiags=0., intDiags=0. }
end

function SpeciesDiagnostics:resetState(tm)
   -- Reset the state indicating if diagnostic has been computed.
   -- This state prevents computing a diagnostic twice when there are dependencies.
   for _, dG in lume.orderedIter(self.diagGroups) do 
      for _, diagNm in ipairs(dG) do self.diags[diagNm].done=false end
   end
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

function SpeciesDiagnostics:writeGridDiagnostics_separateFiles(tm, fr)
   -- Write each grid diagnostic to its own file.
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s_%d.bp", self.name, diagNm, fr), tm, fr, self.writeGhost)
   end
end
function SpeciesDiagnostics:writeIntegratedDiagnostics_separateFiles(tm, fr)
   -- Write each integrated diagnostic to its own file.
   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s.bp", self.name, diagNm), tm, fr)
   end
end

function SpeciesDiagnostics:writeGridDiagnostics_oneFile(tm, fr)
   -- Write all grid diagnostic to a single file.
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      self.gridDiagsToWrite[diagNm] = self.diags[diagNm].field
   end
   self.diagIo:write(self.gridDiagsToWrite, string.format("%s_gridDiagnostics_%d.bp", self.name, fr), tm, fr, self.writeGhost)
end
function SpeciesDiagnostics:writeIntegratedDiagnostics_oneFile(tm, fr)
   -- Write all integrated diagnostic to a single file.
   -- MF 2021/06/01: Not yet available because I think we haven't done this for DynVectors.
   assert(false, "SpeciesDiagnostics:writeIntegratedDiagnostics_oneFile not ready.")
end

function SpeciesDiagnostics:write(tm, fr) self.writeFunc(tm, fr) end

function SpeciesDiagnostics:writeRestartGridDiagnostics_separateFiles(tm, fr)
   -- Write a restart file for each grid diagnostic.
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      diag.field:write(string.format("%s_%s_restart.bp", self.name, diagNm), tm, fr, false)
   end
end
function SpeciesDiagnostics:writeRestartIntegratedDiagnostics_separateFiles(tm, fr)
--   -- Write a restart file for each integrated diagnostic.
--   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
--      local diag = self.diags[diagNm]
--      diag.field:write(string.format("%s_%s_restart.bp", self.name, diagNm), tm, fr, false, false)
--   end
end

function SpeciesDiagnostics:writeRestartGridDiagnostics_oneFile(tm, fr)
   -- Write a single restart file for all grid diagnostics.
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      self.gridDiagsToWrite[diagNm] = self.diags[diagNm].field
   end
   self.diagIo:write(self.gridDiagsToWrite, string.format("%s_gridDiagnostics_restart.bp", self.name), tm, fr, false)
end
function SpeciesDiagnostics:writeRestartIntegratedDiagnostics_oneFile(tm, fr)
   -- Write a single restart file for all integrated diagnostics.
   -- MF 2021/06/01: Not yet available because I think we haven't done this for DynVectors.
   assert(false, "SpeciesDiagnostics:writeRestartIntegratedDiagnostics_oneFile not ready.")
end

function SpeciesDiagnostics:writeRestart(tm, gridFr, intFr) self.writeRestartFunc(tm, gridFr, intFr) end

function SpeciesDiagnostics:readRestartGridDiagnostics_separateFiles()
   local tm, fr
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      local diag = self.diags[diagNm]
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, diagNm))
   end
   return tm, fr
end
function SpeciesDiagnostics:readRestartGridDiagnostics_oneFile()
   for _, diagNm in ipairs(self.diagGroups["grid"]) do 
      self.gridDiagsToRead[diagNm] = self.diags[diagNm].field
   end
   local tm, fr = self.diagIo:read(self.gridDiagsToRead, string.format("%s_gridDiagnostics_restart.bp", self.name))
   return tm, fr
end

function SpeciesDiagnostics:readRestartIntegratedDiagnostics_separateFiles()
   for _, diagNm in ipairs(self.diagGroups["integrated"]) do 
      local diag = self.diags[diagNm]
      diag.field:read(string.format("%s_%s_restart.bp", self.name, diagNm))
   end
end
function SpeciesDiagnostics:readRestartIntegratedDiagnostics_oneFile()
   -- Read from a single restart file for all integrated diagnostics.
   -- MF 2021/06/01: Not yet available because I think we haven't done this for DynVectors.
   assert(false, "SpeciesDiagnostics:readRestartIntegratedDiagnostics_oneFile not ready.")
end

function SpeciesDiagnostics:readRestart()
   local tm, fr = self.readRestartFunc()
   return tm, fr
end

function SpeciesDiagnostics:getDiagTime()
   return self.spentTime.gridDiags + self.spentTime.intDiags
end

return SpeciesDiagnostics
