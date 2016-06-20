-- Gkyl ------------------------------------------------------------------------
--
-- Decomposition algorithms for rectangular grids.
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

-- CartProdDecomp --------------------------------------------------------------
--
-- Decomposition of grid into smaller boxes, specified by number of
-- cuts in each direction
--------------------------------------------------------------------------------

local CartProdDecomp = {}
-- constructor to make new cart-prod decomp
function CartProdDecomp:new(tbl)
   local self = setmetatable({}, CartProdDecomp)
   
   local ndim = #tbl.cuts
   self._cuts = ffi.new("uint16_t[?]", ndim)
   for d = 1, ndim do
      self._cuts[d-1] = tbl.cuts[d]
   end

   self._cellIndex = {} -- start cell indices along each direction
   for d = 1, ndim do
      self._cellIndex[d] = Lin.IntVec(tbl.cuts[d]+1)
   end
   
   return self
end

-- make object callable, and redirect call to the :new method
setmetatable(CartProdDecomp, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
CartProdDecomp.__index = {
   ndim = function (self)
      return #self._cellIndex
   end,
   cuts = function (self, dir)
      return self._cuts[dir-1]
   end,
   decompose = function (self, range) -- decompose range
      for dir = 1, range:ndim() do
	 local cidx = self._cellIndex[dir]
	 cidx[1] = range:lower(dir)

	 local shapes = Lin.IntVec(self:cuts(dir))
	 local baseShape = math.floor(range:shape(dir)/self:cuts(dir))
	 local remCells = range:shape(dir) % self:cuts(dir)	 
	 for c = 1, self:cuts(dir) do
	    shapes[c] = baseShape + (remCells>0 and 1 or 0) -- add extra cell, if any
	    remCells = remCells-1
	 end
	 -- set indices
	 for c = 1, self:cuts(dir) do
	    cidx[c+1] = cidx[c] + shapes[c]
	 end
      end
      return true
   end,
   lower = function (self, dir, c) -- c is cut number
      return self._cellIndex[dir][c]
   end,
   upper = function (self, dir, c) -- c is cut number
      return self._cellIndex[dir][c+1]-1
   end
}

return {
   CartProd = CartProdDecomp,
}
