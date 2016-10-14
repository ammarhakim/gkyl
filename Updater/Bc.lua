-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions specified as BoundaryCondition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Range = require "Lib.Range"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Boundary condition updater
local Bc = {}

function Bc:new (tbl)
   local self = setmetatable({}, Bc)
   Base.setup(self, tbl) -- setup base object
   
   self._grid = assert(tbl.onGrid, "Updater.Bc: Must specify grid to use with 'onGrid''")
   self._dir = assert(tbl.dir, "Updater.Bc: Must specify direction to apply BCs with 'dir'")

   self._edge = assert(
      tbl.edge, "Updater.Bc: Must specify edge to apply BCs with 'edge'. Must be one of 'upper' or 'lower'.")
   if self._edge ~= "lower" and self._edge ~= "upper" then
      error("Updater.Bc: 'edge' must be one of 'lower' or 'upper'. Was " .. self._edge .. " instead")
   end
   self._bcList = assert(
      tbl.boundaryConditions, "Updater.Bc: Must specify boundary conditions to apply with 'boundaryConditions'")

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Bc, { __call = function (self, o) return self.new(self, o) end })

local function getGhostRange(self, global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1 -- for ghost cells on "left"
   else
      lv[self._dir] = global:upper(self._dir)+1 -- for ghost cells on "right"
   end
   return Range.Range(lv, uv)
end

local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir, edge = self._dir, self._edge
   local global, globalExt = qOut:globalRange(), qOut:globalExtRange()
   local localExtRange = qOut:localExtRange()
   local ghost = localExtRange:intersect(
      getGhostRange(self, global, globalExt)) -- range spanning ghost cells
   
   local qG, qS = qOut:get(0), qOut:get(0) -- get pointers to (re)use inside inner loop [G: Ghost, S: Skin]
   
   local indexer = qOut:genIndexer()
   for idx in ghost:colMajorIter() do -- loop, applying BCs
      local idxs = idx:copy()
      idxs[dir] = self._edge == "lower" and global:lower(dir)  or global:upper(dir)
      qOut:fill(indexer(idx), qG); qOut:fill(indexer(idxs), qS)
      for _, bc in ipairs(self._bcList) do -- loop over each BC
	 bc(dir, tCurr+dt, nil, qS, qG) -- TODO: PASS COORDINATES
      end
   end
   return true, 1e37
end

-- Methods for BC object
Bc.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   Bc = Bc
}
