-- Gkyl ------------------------------------------------------------------------
--
-- Modal Maxium-Order elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- CartModalMaxOrder -----------------------------------------------------------
--
-- Modal maximal-order basis set
--------------------------------------------------------------------------------

local CartModalMaxOrder = {}
function CartModalMaxOrder:new(tbl)
   local self = setmetatable({}, CartModalMaxOrder)

   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalMaxOrder: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalMaxOrder: Must specify polynonial order with 'polyOrder'")

   if (self._polyOrder < 0) or (self._polyOrder > 4) then
      assert(false, "Polynomial order must be between 0 and 4")
   end

   self._numBasis = 1

   -- number of basis is = (p+d)! / p! d! 
   if self._polyOrder > 0 then
      if (self._ndim == 1) then
	 self._numBasis = self._polyOrder+1
      elseif (self._ndim == 2) then
	 self._numBasis = (self._polyOrder+2)*(self._polyOrder+1)/2
      elseif (self._ndim == 3) then
	 self._numBasis = (self._polyOrder+3)*(self._polyOrder+2)*(self._polyOrder+1)/6
      elseif (self._ndim == 4) then
	 self._numBasis = (self._polyOrder+4)*(self._polyOrder+3)*(self._polyOrder+2)*(self._polyOrder+1)/24
      elseif (self._ndim == 5) then
	 self._numBasis = (self._polyOrder+5)*(self._polyOrder+4)*(self._polyOrder+3)*(self._polyOrder+2)*(self._polyOrder+1)/120
      end
   end

   local _m = nil -- to store module with evaluation code
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalBasis1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalMaxOrderBasis2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalMaxOrderBasis3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalMaxOrderBasis4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalMaxOrderBasis5d"
   end

   self._evalBasisFunc = _m[self._polyOrder]
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartModalMaxOrder, { __call = function (self, o) return self.new(self, o) end })

CartModalMaxOrder.__index = {
   ndim = function (self)
      return self._ndim
   end,
   polyOrder = function (self)
      return self._polyOrder
   end,
   numBasis = function (self)
      return self._numBasis
   end,
   evalBasis = function (self, z, b)
      return self._evalBasisFunc(z, b)
   end,
}

return {
   CartModalMaxOrder = CartModalMaxOrder   
}
