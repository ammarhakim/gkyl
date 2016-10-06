-- Gkyl ------------------------------------------------------------------------
--
-- Compute and store neighbor information for communication in parallel
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"

local CartDecompNeigh = {}
-- constructor to make new neighbor object
function CartDecompNeigh:new(decomp)
   local self = setmetatable({}, CartDecompNeigh)
   self._decomp = decomp -- decomposition object
   self._neighData = {} -- neighbor data
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartDecompNeigh, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
CartDecompNeigh.__index = {
   calcFaceCommNeigh = function (self, lowerGhost, upperGhost)
      self._neighData = {} -- clear out existing data
      local ndim = self._decomp:ndim()
      local numSubDomains = self._decomp:numSubDomains()
      for kme = 1, numSubDomains do
	 local nlst = {} -- List of neighbors
	 for d = 1, ndim  do
	    -- expand sub-domain in direction `d`
	    local expSubDom = self._decomp:subDomain(kme):extendDir(d, lowerGhost, upperGhost)
	    -- loop over all other sub-domains and intersect
	    for ku = 1, numSubDomains do
	       if ku == kme then goto continue end -- no self-intersections
	       if not expSubDom:isIntersectionEmpty(self._decomp:subDomain(ku)) then
		  table.insert(nlst, ku) -- insert subDomain index into list of neighbors
	        end
	       ::continue::
	    end
	 end
	 self._neighData[kme] = nlst
      end
   end,
   calcAllCommNeigh = function (self, lowerGhost, upperGhost)
      self._neighData = {} -- clear out existing data
      local ndim = self._decomp:ndim()
      local numSubDomains = self._decomp:numSubDomains()
      for kme = 1, numSubDomains do
	 local nlst = {} -- List of neighbors
	 -- expand sub-domain
	 local expSubDom = self._decomp:subDomain(kme):extend(lowerGhost, upperGhost)
	 -- loop over all other sub-domains and intersect
	 for ku = 1, numSubDomains do
	    if ku == kme then goto continue end -- no self-intersections
	    if not expSubDom:isIntersectionEmpty(self._decomp:subDomain(ku)) then
	       table.insert(nlst, ku) -- insert subDomain index into list of neighbors
	    end
	    ::continue::
	 end
	 self._neighData[kme] = nlst
      end
   end,
   neighborData = function(self, k)
      return self._neighData[k] and self._neighData[k] or {}
   end
}

return CartDecompNeigh
