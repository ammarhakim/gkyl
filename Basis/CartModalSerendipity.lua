-- Gkyl ------------------------------------------------------------------------
--
-- Modal Serendipity elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- CartModalSerendipity -----------------------------------------------------------
--
-- Modal serendipity basis set. See: Found. Comput. Math 2011 11:337 Arnold & Awanou
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local function numBasis(ndim, polyOrder)
   local numBasis_d2 = {4, 8, 12, 17}
   local numBasis_d3 = {8, 20, 32, 50}
   local numBasis_d4 = {16, 48, 80, 136}
   local numBasis_d5 = {32, 112, 192, 352}

   local nbasis = 1
   if polyOrder > 0 then
      if (ndim == 1) then
	 nbasis = polyOrder+1
      elseif (ndim == 2) then
	 nbasis = numBasis_d2[polyOrder]
      elseif (ndim == 3) then
	 nbasis = numBasis_d3[polyOrder]
      elseif (ndim == 4) then
	 nbasis = numBasis_d4[polyOrder]
      elseif (ndim == 5) then
	 nbasis = numBasis_d5[polyOrder]
      end
   end
   return nbasis
end

local CartModalSerendipity = {}
function CartModalSerendipity:new(tbl)
   local self = setmetatable({}, CartModalSerendipity)

   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalSerendipity: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalSerendipity: Must specify polynonial order with 'polyOrder'")

   if (self._polyOrder < 0) or (self._polyOrder > 4) then
      assert(false, "Polynomial order must be between 0 and 4")
   end

   self._numBasis = numBasis(self._ndim, self._polyOrder)
   self._numSurfBasis = numBasis(self._ndim-1, self._polyOrder)

   local _m = nil -- to store module with evaluation code
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalBasis1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalSerendipBasis2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalSerendipBasis3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalSerendipBasis4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalSerendipBasis5d"
   end

   self._evalBasisFunc = _m[self._polyOrder]
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartModalSerendipity, { __call = function (self, o) return self.new(self, o) end })

CartModalSerendipity.__index = {
   id = function (self)
      return "serendipity"
   end,   
   ndim = function (self)
      return self._ndim
   end,
   polyOrder = function (self)
      return self._polyOrder
   end,
   numBasis = function (self)
      return self._numBasis
   end,
   numSurfBasis = function (self)
      return self._numSurfBasis
   end,   
   evalBasis = function (self, z, b)
      return self._evalBasisFunc(z, b)
   end,
}

return {
   CartModalSerendipity = CartModalSerendipity   
}
