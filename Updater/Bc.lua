-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions specified as BoundaryCondition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"

-- Boundary condition updater.
local Bc = Proto(UpdaterBase)

function Bc:init(tbl)
   Bc.super.init(self, tbl) -- Setup base object.

   self._isFirst = true -- Will be reset first time _advance() is called.

   self._grid = assert(tbl.onGrid, "Updater.Bc: Must specify grid to use with 'onGrid''")
   self._dir  = assert(tbl.dir, "Updater.Bc: Must specify direction to apply BCs with 'dir'")

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

   self._ghostRangeDecomp = nil -- Will be constructed on first call to advance.

   self.hasExtFld = xsys.pickBool(tbl.hasExtFld, false)
end

function Bc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1 -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1 -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function Bc:_advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")
   local qIn  = inFld[1]

   local dir, edge = self._dir, self._edge
   local vdir      = self._vdir
   local global    = qOut:globalRange()

   if self._isFirst then
      local globalExt     = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      local ghostRng      = localExtRange:intersect(
   	 self:getGhostRange(global, globalExt)) -- Range spanning ghost cells.
      if self._skinLoop == "integrate" then
   	 self._skin = qOut:localRange():selectLast(self._vdim)
      end
      -- Decompose ghost region into threads.
      self._ghostRangeDecomp = LinearDecomp.LinearDecompRange {
   	 range = ghostRng, numSplit = grid:numSharedProcs() }
   end

   local qG, qS = qOut:get(1), qOut:get(1) -- Get pointers to (re)use inside inner loop [G: Ghost, S: Skin].
   if self.hasExtFld then qS = qIn:get(1) end
   local idxS = Lin.IntVec(grid:ndim()) -- Prealloc this.
   local indexer = qOut:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   for idxG in self._ghostRangeDecomp:rowMajorIter(tId) do -- Loop, applying BCs.
      idxG:copyInto(idxS)
      -- If an in-field is specified the same indexes are used (gS
      -- points to the ghost layer of the in-field); otherwise, move
      -- the ghost index to point into the skin layer.
      if not self.hasExtFld then
	 idxS[dir] = edge == "lower" and global:lower(dir) or global:upper(dir)
	 if self._skinLoop == "flip" then
	    idxS[vdir] = global:upper(vdir) + 1 - idxS[vdir]
	 end
      end
      qOut:fill(indexer(idxG), qG) 
      if self._skinLoop == "integrate" then 
   	 for c = 1, self._numComponents do qG[c] = 0 end

         for idx in self._skin:rowMajorIter() do
   	    for d = 1, self._vdim do idxS[self._cdim + d] = idx[d] end
   	    qOut:fill(indexer(idxS), qS)
            for _, bc in ipairs(self._bcList[1]) do
               bc(dir, tCurr, idxS, qS, qG, self._bcList[2])
            end
         end
      else
	 if not self.hasExtFld then
	    qOut:fill(indexer(idxS), qS)
	 else
	    qIn:fill(indexer(idxS), qS)
	 end
         for _, bc in ipairs(self._bcList[1]) do
            bc(dir, tCurr, idxS, qS, qG, self._bcList[2]) -- TODO: PASS COORDINATES.
         end
      end
   end

   self._isFirst = false
end

function Bc:getDir()
   return self._dir
end

return Bc
