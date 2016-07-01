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

function _M.Copy(tbl)
   local cIdx = getComponents(tbl) -- components to apply to
   local n = #cIdx -- number of components
   return function (dir, tm, xc, qin, qbc)
      for i = 1, n do
	 qbc[cIdx[i]] = qin[cIdx[i]]
      end
   end
end

return _M


