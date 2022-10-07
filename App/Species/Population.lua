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

--function Population:beginCartFieldTransfer(srcDestRank, fldIn, fldOut)
--   local tag = 32
--   local recvReq = Mpi.Irecv(fldOut:dataPointer(), fldOut:size(), fldOut:elemCommType(), srcDestRank, tag, self.comm)
--   local sendReq = Mpi.Isend(fldIn:dataPointer(), fldIn:size(), fldIn:elemCommType(), srcDestRank, tag, self.comm)
--   return recvReq, sendReq
--end
--
--function Population:initCrossSpeciesCoupling()
--   -- Create a double nested table of colliding species.
--   -- In this table we will encode information about that collition such as:
--   --   * does the collision take place?
--   -- Other features of a collision may be added in the future.
--   self.collPairs = {}
--   for sN, _ in lume.orderedIter(species) do
--      self.collPairs[sN] = {}
--      for sO, _ in lume.orderedIter(species) do
--         self.collPairs[sN][sO] = {}
--         -- Need next below because species[].collisions is created as an empty table.
--         if species[sN].collisions and next(species[sN].collisions) then
--            for collNm, _ in pairs(species[sN].collisions) do
--               -- This species collides with someone.
--               self.collPairs[sN][sO].on = lume.any(species[sN].collisions[collNm].collidingSpecies,
--                                                    function(e) return e==sO end)
--            end
--         else
--            -- This species does not collide with anyone.
--            self.collPairs[sN][sO].on = false
--         end
--      end
--   end
--
----   -- Create list of ranks we need to send/recv local threeMoments to/from.
----   self.threeMomentsXfer = {}
----   self.threeMomentsXfer.to, self.threeMoments.from = {}, {}
----   self.threeMomentsXfer.sendReq, self.threeMoments.recvReq = {}, {}
----   for sN, _ in self.iterLocal() do
----      for sO, v in pairs(self.collPairs[sN]) do
----         local sOrank = self:getSpeciesOwner(sO)
----         if sO~=sN and v.on and (not self:isSpeciesMine(sO)) then
----            if (not lume.any(self.threeMomentsXfer.to, function(e) return e==sOrank end)) then 
----               self.threeMomentsXfer.to[sN] = sOrank 
----               self.threeMomentsXfer.sendReq[sN] = Mpi.Request()
----            end
----            self.threeMomentsXfer.from[sO] = sOrank 
----            self.threeMomentsXfer.recvReq[sO] = Mpi.Request()
----         end
----      end
----   end
----
----   -- Allocate space threeMoments in species we collide with but not in this rank.
----   for sLnm, _ in self.iterLocal() do
----      for sOnm, v in pairs(self.collPairs[sLnm]) do
----         local sO = self.species[sOnm]
----         if v.on and (not self:isSpeciesMine(sOnm)) and sO.threeMoments==nil then
----            sO.threeMoments = sO:allocVectorMoment(sO.udim+2)
----         end
----      end
----   end
--end
--
--function Population:beginCouplingMomentsXfer()
--   -- When parallelizing over species, and if we don't own other species
--   -- this species collides with, we'll need to communicate threeMoments.
--
--   -- Post receives for every species not in this rank.
--   local tag = 32
--   for nm, rank in pairs(self.threeMomentsXfer.from) do
--      local recvFld = self.species[nm].threeMoments
--      local req = self.threeMomentsXfer.recvReq[nm]
--      local err = Mpi.Irecv(recvFld:dataPointer(), recvFld:size(),
--                            recvFld:elemCommType(), rank, tag, self.comm, req)
--   end
--   -- Post sends for every species in this rank.
--   for nm, rank in pairs(self.threeMomentsXfer.to) do
--      local sendFld = self.species[nm].threeMoments
--      local req = self.threeMomentsXfer.sendReq[nm]
--      local err = Mpi.Isend(sendFld:dataPointer(), sendFld:size(),
--                            sendFld:elemCommType(), rank, tag, self.comm, req)
--   end
--end
--
--function Population:beginSpeciesXfer(speciesNm, fld)
--   -- Begin nonblocking communicaiton of field 'fld' from species named 'speciesNm'. 
--
--   -- Post receives for every species not in this rank.
--   local tag = 32
--   for nm, rank in pairs(self.threeMomentsXfer.from) do
--      local recvFld = self.species[nm].threeMoments
--      local req = self.threeMomentsXfer.recvReq[nm]
--      local err = Mpi.Irecv(recvFld:dataPointer(), recvFld:size(),
--                            recvFld:elemCommType(), rank, tag, self.comm, req)
--   end
--   -- Post sends for every species in this rank.
--   for nm, rank in pairs(self.threeMomentsXfer.to) do
--      local sendFld = self.species[nm].threeMoments
--      local req = self.threeMomentsXfer.sendReq[nm]
--      local err = Mpi.Isend(sendFld:dataPointer(), sendFld:size(),
--                            sendFld:elemCommType(), rank, tag, self.comm, req)
--   end
--end

return Population
