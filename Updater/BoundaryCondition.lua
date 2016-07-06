-- Gkyl ------------------------------------------------------------------------
--
-- Boundary condition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local _M = {}

-- Helper function to extract list of components
local function getComponents(tbl)
   local c = assert(tbl.components, "BoundaryCondition: Must specify components with 'components' parameter")
   local cIdx = Lin.IntVec(#c)
   for i = 1, #c do
      cIdx[i] = c[i]
   end
   return cIdx
end

-- copy BCs
function _M.Copy(tbl)
   local cIdx = getComponents(tbl) -- components to apply to
   local n = #cIdx -- number of components
   return function (dir, tm, xc, qin, qbc)
      for i = 1, n do
	 qbc[cIdx[i]] = qin[cIdx[i]]
      end
   end
end

-- Resuffle indices for various direction normal vectors. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int32_t[4]", 0, 1, 2, 3),
   new("int32_t[4]", 0, 2, 3, 1),
   new("int32_t[4]", 0, 3, 1, 2)
}

-- zero normal
function _M.ZeroNormal(tbl)
   local cIdx = getComponents(tbl) -- components to apply to
   local n = #cIdx -- number of components
   assert(n == 3, "ZeroNormal BC can only be applied to vector of length 3")
   return function (dir, tm, xc, qin, qbc)
      local d = dirShuffle[dir]
      qbc[cIdx[d[1]]] = -qin[cIdx[d[1]]]
      qbc[cIdx[d[2]]] = qin[cIdx[d[2]]]
      qbc[cIdx[d[3]]] = qin[cIdx[d[3]]]
   end
end

-- zero tangent
function _M.ZeroTangent(tbl)
   local cIdx = getComponents(tbl) -- components to apply to
   local n = #cIdx -- number of components
   assert(n == 3, "ZeroTangent BC can only be applied to vector of length 3")
   return function (dir, tm, xc, qin, qbc)
      local d = dirShuffle[dir]
      qbc[cIdx[d[1]]] = qin[cIdx[d[1]]]
      qbc[cIdx[d[2]]] = -qin[cIdx[d[2]]]
      qbc[cIdx[d[3]]] = -qin[cIdx[d[3]]]
   end
end

return _M


