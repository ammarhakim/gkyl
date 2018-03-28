-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions specified as BoundaryCondition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Lin = require "Lib.Linalg"
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

   self._skinLoop = tbl.skinLoop
   if self._skinLoop == nil then
      self._skinLoop = "pointwise"
   end
   if self._skinLoop ~= "pointwise" and self._skinLoop ~= "flip" and self._skinLoop ~= "integrate" then
      error("Updater.Bc: 'skinLoop' must one of 'pointwise', 'flip', or 'integrate'."
	       .. " Was '" .. self._skinLoop .. "' instead")
   end
   if self._skinLoop == "flip" or self._skinLoop == "integrate" then
      self._cdim = assert(
	 tbl.cdim,  "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
   end

   self._ghost = nil -- store ghost cells
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

function Bc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir, edge = self._dir, self._edge
   local global = qOut:globalRange()

   if self._isFirst then
      -- compute ghost cells first time around only
      local globalExt = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      self._ghost = localExtRange:intersect(
	 self:getGhostRange(global, globalExt)) -- range spanning ghost cells
   end

   local qG, qS = qOut:get(1), qOut:get(1) -- get pointers to (re)use inside inner loop [G: Ghost, S: Skin]
   local idxS = Lin.IntVec(grid:ndim()) -- prealloc this
   local indexer = qOut:genIndexer()

   if self._skinLoop == "pointwise" then
      for idx in self._ghost:colMajorIter() do -- loop, applying BCs
	 idx:copyInto(idxS)
	 idxS[dir] = edge == "lower" and global:lower(dir)  or global:upper(dir)
	 qOut:fill(indexer(idx), qG); qOut:fill(indexer(idxS), qS)
	 for _, bc in ipairs(self._bcList) do -- loop over each BC
	    bc(dir, tCurr+dt, nil, qS, qG) -- TODO: PASS COORDINATES
	 end
      end
   elseif self._skinLoop == "flip" then
      local vdir = dir + self._cdim
      for idx in self._ghost:colMajorIter() do -- loop, applying BCs
	 idx:copyInto(idxS)
	 idxS[dir] = edge == "lower" and global:lower(dir)  or global:upper(dir)
	 idxS[vdir] = global:upper(vdir) + 1 - idxS[vdir]
	 qOut:fill(indexer(idx), qG); qOut:fill(indexer(idxS), qS)
	 for _, bc in ipairs(self._bcList) do -- loop over each BC
	    bc(dir, tCurr+dt, nil, qS, qG) -- TODO: PASS COORDINATES
	 end
      end
   end

   self._isFirst = false
   return true, GKYL_MAX_DOUBLE
end

return Bc
