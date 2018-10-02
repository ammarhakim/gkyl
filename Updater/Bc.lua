-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions specified as BoundaryCondition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"

-- Boundary condition updater
local Bc = Proto(UpdaterBase)

function Bc:init(tbl)
   Bc.super.init(self, tbl) -- setup base object

   self._isFirst = true -- will be reset first time _advance() is called

   self._grid = assert(tbl.onGrid, "Updater.Bc: Must specify grid to use with 'onGrid''")
   self._dir = assert(tbl.dir, "Updater.Bc: Must specify direction to apply BCs with 'dir'")

   self._edge = assert(
      tbl.edge, "Updater.Bc: Must specify edge to apply BCs with 'edge'. Must be one of 'upper', or 'lower'.")
   if self._edge ~= "lower" and self._edge ~= "upper" then
      error("Updater.Bc: 'edge' must be one of 'lower' or 'upper'. Was " .. self._edge .. " instead")
   end
   self._bcList = assert(
      tbl.boundaryConditions, "Updater.Bc: Must specify boundary conditions to apply with 'boundaryConditions'")

   self._skinLoop = tbl.skinLoop and tbl.skinLoop or "pointwise"
   if self._skinLoop == "flip" or self._skinLoop == "integrate" then
      self._cdim = assert(
	 tbl.cdim,
	 "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
   end
   if self._skinLoop == "flip" then
      self._vdir = assert(
	 tbl.vdir,
	 "Updater.Bc: Must specify velocity direction to flip with 'vdir'")
   end
   if self._skinLoop == "integrate" then
      self._vdim = assert(
	 tbl.vdim,
	 "Updater.Bc: Must specify velocity space dimensions to apply with 'vdim'")
      self._numComps = assert(
	 tbl.numComps,
	 "Updater.Bc: Must specify the number of components of the field using 'numComps'")
   end

   self._ghost = nil -- store ghost cells
   self._ghostRangeDecomp = nil -- ghost cell range shared decomposition
   self._localStartIdx, self._localNumBump = nil, nil -- local start index and bump for shared decomp of ghost range

   -- get the index of the shared processor
   self._shmIndex = Mpi.Comm_rank(self:getShmComm())+1 -- our local index on SHM comm (one more than rank)
end

function Bc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1 -- for ghost cells on "left"
   else
      lv[self._dir] = global:upper(self._dir)+1 -- for ghost cells on "right"
   end
   return Range.Range(lv, uv)
end

function Bc:getVelocityRange(localRange)
   local lv, uv = {}, {}
   for d = 1, self._vdim do
      lv[d]= localRange:lower(self._cdim + d)
      uv[d]= localRange:upper(self._cdim + d)
   end
   return Range.Range(lv, uv)
end

function Bc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir, edge = self._dir, self._edge
   local vdir = self._vdir
   local global = qOut:globalRange()

   -- get ahold of communicators
   local shmComm = self:getShmComm()
   local worldComm = self:getWorldComm()

   if self._isFirst then
      -- compute ghost cells first time around only
      local globalExt = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      self._ghost = localExtRange:intersect(
	 self:getGhostRange(global, globalExt)) -- range spanning ghost cells
      self._ghostRangeDecomp = LinearDecomp.LinearDecompRange {
         range = self._ghost,
         numSplit = Mpi.Comm_size(shmComm)
      }
      -- store start index and size handled by local SHM-rank for local range
      self._localStartIdx, self._localNumBump = self._ghostRangeDecomp:colStartIndex(self._shmIndex), self._ghostRangeDecomp:shape(self._shmIndex)
      if self._skinLoop == "integrate" then
	 self._skin = self:getVelocityRange(qOut:localRange())
      end
   end

   local qG, qS = qOut:get(1), qOut:get(1) -- get pointers to (re)use inside inner loop [G: Ghost, S: Skin]
   local idxS = Lin.IntVec(grid:ndim()) -- prealloc this
   local indexer = qOut:genIndexer()

   -- barrier before shared loop
   Mpi.Barrier(worldComm)

   for idxG in self._ghost:colMajorIter(self._localStartIdx, self._localNumBump) do -- loop, applying BCs
      idxG:copyInto(idxS)
      idxS[dir] = edge == "lower" and global:lower(dir) or global:upper(dir)
      if self._skinLoop == "flip" then
	 idxS[vdir] = global:upper(vdir) + 1 - idxS[vdir]
      end
      qOut:fill(indexer(idxG), qG) 
      if self._skinLoop == "integrate" then 
	 -- clear the ghost cells before accumulating
	 for c = 1, self._numComponents do qG[c] = 0 end

         for idx in self._skin:colMajorIter() do
	    for d = 1, self._vdim do idxS[self._cdim + d] = idx[d] end
	    qOut:fill(indexer(idxS), qS)
            for _, bc in ipairs(self._bcList) do -- loop over each BC
               bc(dir, tCurr+dt, idxS, qS, qG)
            end
         end
      else
         qOut:fill(indexer(idxS), qS)
         for _, bc in ipairs(self._bcList) do -- loop over each BC
            bc(dir, tCurr+dt, idxS, qS, qG) -- TODO: PASS COORDINATES
         end
      end
   end

   self._isFirst = false
   return true, GKYL_MAX_DOUBLE
end

return Bc
