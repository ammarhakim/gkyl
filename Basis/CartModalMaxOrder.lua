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

   self._numBasis = numBasis(self._ndim, self._polyOrder)
   self._numSurfBasis = numBasis(self._ndim-1, self._polyOrder)

   local _m = nil -- to store module with evaluation code
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      assert(self._polyOrder <= 7, "For 1D polynomial order must be either 1, 2, 3, 4, 5, 6, or 7")
      _m = require "Basis._data.ModalBasis1d"
   elseif (self._ndim == 2) then
      assert(self._polyOrder <= 7, "For 2D polynomial order must be either 1, 2, 3, 4, 5, 6, or 7")
      _m = require "Basis._data.ModalMaxOrderBasis2d"
   elseif (self._ndim == 3) then
      assert(self._polyOrder <= 6, "For 3D polynomial order must be either 1, 2, 3, 4, 5, or 6")
      _m = require "Basis._data.ModalMaxOrderBasis3d"
   elseif (self._ndim == 4) then
      assert(self._polyOrder <= 5, "For 4D polynomial order must be either 1, 2, 3, 4, or 5")
      _m = require "Basis._data.ModalMaxOrderBasis4d"
   elseif (self._ndim == 5) then
      assert(self._polyOrder <= 4, "For 5D polynomial order must be either 1, 2, 3, or 4")
      _m = require "Basis._data.ModalMaxOrderBasis5d"
   elseif (self._ndim == 6) then
      assert(self._polyOrder <= 3, "For 6D polynomial order must be either 1, 2, or 3")
      _m = require "Basis._data.ModalMaxOrderBasis6d"
   end

   self._evalBasisFunc = _m[self._polyOrder] -- function to evaluate basis functions

   _m = nil -- to store module with flip-sign method
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalBasisFlipSign1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalMaxOrderBasisFlipSign2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalMaxOrderBasisFlipSign3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalMaxOrderBasisFlipSign4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalMaxOrderBasisFlipSign5d"
   elseif (self._ndim == 6) then
      _m = require "Basis._data.ModalMaxOrderBasisFlipSign6d"
   end

   self._flipSign = _m[self._polyOrder] -- function to flip sign
   
end

function CartModalMaxOrder:id() return "maximal-order" end
function CartModalMaxOrder:ndim() return self._ndim end
function CartModalMaxOrder:polyOrder() return self._polyOrder end
function CartModalMaxOrder:numBasis() return self._numBasis end
function CartModalMaxOrder:numSurfBasis() return self._numSurfBasis end
function CartModalMaxOrder:evalBasis(z, b) return self._evalBasisFunc(z, b) end
function CartModalMaxOrder:flipSign(dir, fIn, fOut) self._flipSign(dir, fIn, fOut) end

return {
   CartModalMaxOrder = CartModalMaxOrder   
}
