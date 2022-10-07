-- Gkyl ------------------------------------------------------------------------
--
-- Population object (a collection of species).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Mpi   = require "Comm.Mpi"
local Lin   = require "Lib.Linalg"
local lume  = require "Lib.lume"

-- Empty shell species base class.
local Population = Proto()

-- Functions that must be defined by subclasses.
function Population:init(tbl)
   self.comm = assert(tbl.comm, "App.Population: must specify the species communicator with 'comm'.")
   self.rank = Mpi.Comm_rank(self.comm)

   self.decompCuts = Mpi.Comm_size(self.comm)

   self.species = {}
end

function Population:decompose()
   -- Distribute species across the species communicator. IMPORTANT: call this after
   -- the species names have been recorded and lume.setOrder(species) has been called.
   local numSpecies = 0
   for nm, _ in pairs(self.species) do numSpecies=numSpecies+1 end

   local shapes     = Lin.IntVec(self.decompCuts)
   local baseShape  = math.floor(numSpecies/self.decompCuts)

   local remSpec = numSpecies % self.decompCuts -- Extra species left over.
   for c = 1, self.decompCuts do
      shapes[c] = baseShape + (remSpec>0 and 1 or 0) -- Add extra species if any remain.
      remSpec = remSpec-1
   end
   local starts = Lin.IntVec(self.decompCuts)
   starts[1] = 1
   for c = 2, self.decompCuts do starts[c] = starts[c-1] + shapes[c-1] end

   -- Create a table of the species in this rank.
   self.mySpeciesNames = {}
   local keys = getmetatable(self.species)
   for i = starts[self.rank+1], starts[self.rank+1]+shapes[self.rank+1]-1 do
      table.insert(self.mySpeciesNames, keys[i])
   end

   -- Set local and global iterators in population table.
   self.iterLocal = function()
      local idx, keys = 0, getmetatable(self.species)
      local keysLocal = self.mySpeciesNames
      local n = table.getn(keysLocal)
      return function()
         idx = idx+1
         if idx <= n then return keysLocal[idx], self.species[keysLocal[idx]] end
      end
   end

   self.iterGlobal = function() return lume.orderedIter(self.species) end

   -- So we don't have to search into mySpeciesName every time,
   -- create a key index table indicating if species is in this rank.
   self.isMySpecies = {}
   for nm, _ in self.iterGlobal() do
      self.isMySpecies[nm] = lume.any(self.mySpeciesNames, function(e) return e==nm end)
   end

   -- Create a table of ranks containing other species:
   self.speciesOwnerRank = {}
   for r = 1, self.decompCuts do
      for i = starts[r], starts[r]+shapes[r]-1 do self.speciesOwnerRank[keys[i]] = r-1 end
   end

   -- Create an iterator over the species not in this rank.
   local foreignSpeciesNames = {}
   for nm, _ in self.iterGlobal() do
      if not self.isMySpecies[nm] then table.insert(foreignSpeciesNames, nm) end
   end
   self.iterForeign = function()
      local idx, keysForeign = 0, foreignSpeciesNames
      local n = table.getn(keysForeign)
      return function()
         idx = idx+1
         if idx <= n then return keysForeign[idx], self.species[keysForeign[idx]] end
      end
   end
end

function Population:isSpeciesMine(speciesName)
   return self.isMySpecies[speciesName]
end

function Population:getMass(speciesName)
   return self.species[speciesName]:getMass()
end

function Population:getCharge(speciesName)
   return self.species[speciesName]:getCharge()
end

function Population:getSpeciesOwner(speciesName)
   return self.speciesOwnerRank[speciesName]
end

function Population:getComm() return self.comm end
function Population:getSpecies() return self.species end

function Population:speciesXferField_begin(xferObj, fld, tag)
   -- Initiate the communication of a CartField between species. This is to be called
   -- by a species or collision App, which has to specify a transfer object containing
   -- the destination and source ranks, and the necessary MPI requests/statuses.

   -- Post receive from rank owning this species.
   for _, rank in ipairs(xferObj.srcRank) do
      local _ = Mpi.Irecv(fld:dataPointer(), fld:size(), fld:elemCommType(),
                          rank, tag, self:getComm(), xferObj.recvReqStat)
   end
   -- Post sends to species that don't have this species.
   for _, rank in ipairs(xferObj.destRank) do
      local _ = Mpi.Isend(fld:dataPointer(), fld:size(), fld:elemCommType(),
                          rank, tag, self:getComm(), xferObj.sendReqStat)
   end
end

function Population:speciesXferField_waitRecv(xferObj)
   for _, rank in ipairs(xferObj.srcRank) do
      local _ = Mpi.Wait(xferObj.recvReqStat, xferObj.recvReqStat)
   end
end

function Population:speciesXferField_waitSend(xferObj)
   for _, rank in ipairs(xferObj.destRank) do
      local _ = Mpi.Wait(xferObj.sendReqStat, xferObj.sendReqStat)
   end
end

return Population
