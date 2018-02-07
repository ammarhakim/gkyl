-- Gkyl ------------------------------------------------------------------------
--
-- Modal Maximum-Order elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"

-- CartModalMaxOrder -----------------------------------------------------------
--
-- Modal maximal-order basis set
--------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local function numBasis(ndim, polyOrder)
   local nbasis = 1
   -- number of basis is = (p+d)! / p! d! 
   if (ndim == 0) then
      nbasis = 1
   elseif (ndim == 1) then
      nbasis = polyOrder+1
   elseif (ndim == 2) then
      nbasis = (polyOrder+2)*(polyOrder+1)/2
   elseif (ndim == 3) then
      nbasis = (polyOrder+3)*(polyOrder+2)*(polyOrder+1)/6
   elseif (ndim == 4) then
      nbasis = (polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/24
   elseif (ndim == 5) then
      nbasis = (polyOrder+5)*(polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/120
   elseif (ndim == 6) then
      nbasis = (polyOrder+6)*(polyOrder+5)*(polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/720
   end
   return nbasis
end

local CartModalMaxOrder = Proto()
function CartModalMaxOrder:init(tbl)
   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalMaxOrder: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalMaxOrder: Must specify polynonial order with 'polyOrder'")

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
      _m = require "Basis._data.ModalMaxOrderBasis2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalMaxOrderBasis3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalMaxOrderBasis4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalMaxOrderBasis5d"
   elseif (self._ndim == 6) then
      assert(self._polyOrder <= 2, "For 6D polynomial order must be either 1 or 2")
      _m = require "Basis._data.ModalMaxOrderBasis6d"
   end

   self._evalBasisFunc = _m[self._polyOrder] -- function to evaluate basis functions
end

function CartModalMaxOrder:id() return "maximal-order" end
function CartModalMaxOrder:ndim() return self._ndim end
function CartModalMaxOrder:polyOrder() return self._polyOrder end
function CartModalMaxOrder:numBasis() return self._numBasis end
function CartModalMaxOrder:numSurfBasis() return self._numSurfBasis end
function CartModalMaxOrder:evalBasis(z, b) return self._evalBasisFunc(z, b) end

return {
   CartModalMaxOrder = CartModalMaxOrder   
}
