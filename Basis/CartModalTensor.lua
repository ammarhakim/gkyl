-- Gkyl ------------------------------------------------------------------------
--
-- Modal Tensor elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"

-- CartModalTensor -----------------------------------------------------------
--
-- Modal Tensor basis set.
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local function numBasis(ndim, polyOrder)
   local nbasis = 1
   -- number of basis is = (p+1)^d
   if (ndim == 0) then
      nbasis = 1
   elseif (ndim == 1) then
      nbasis = (polyOrder+1)
   elseif (ndim == 2) then
      nbasis = (polyOrder+1)*(polyOrder+1)
   elseif (ndim == 3) then
      nbasis = (polyOrder+1)*(polyOrder+1)*(polyOrder+1)
   elseif (ndim == 4) then
      nbasis = (polyOrder+1)*(polyOrder+1)*(polyOrder+1)*(polyOrder+1)
   end
   return nbasis
end

local CartModalTensor = Proto()
function CartModalTensor:init(tbl)
   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalTensor: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalTensor: Must specify polynonial order with 'polyOrder'")

   self._numBasis = numBasis(self._ndim, self._polyOrder)
   self._numSurfBasis = numBasis(self._ndim-1, self._polyOrder)

   local _m = nil -- to store module with evaluation code
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      assert(self._polyOrder <= 7, "For 1D polynomial order must be either 1, 2, 3, 4, 5, 6, or 7")
      _m = require "Basis._data.ModalBasis1d"
   elseif (self._ndim == 2) then
      assert(self._polyOrder <= 4, "For 2D polynomial order must be either 1, 2, 3, or 4")
      _m = require "Basis._data.ModalTensorBasis2d"
   elseif (self._ndim == 3) then
      assert(self._polyOrder <= 3, "For 3D polynomial order must be either 1, 2, or 3")
      _m = require "Basis._data.ModalTensorBasis3d"
   elseif (self._ndim == 4) then
      assert(self._polyOrder <= 2, "For 4D polynomial order must be either 1 or 2")
      _m = require "Basis._data.ModalTensorBasis4d"
   elseif (self._ndim == 5) then
   -- piecewise linear Serendipity and piecewise linear Tensor are equivalent
      assert(self._polyOrder <= 1, "For 5D polynomial order must be 1")
      _m = require "Basis._data.ModalSerendipBasis5d"
   elseif (self._ndim == 6) then
   -- piecewise linear Serendipity and piecewise linear Tensor are equivalent
      assert(self._polyOrder <= 1, "For 6D polynomial order must be 1")
      _m = require "Basis._data.ModalSerendipBasis6d"

   end

   self._evalBasisFunc = _m[self._polyOrder] -- function to evaluate basis functions

   _m = nil -- to store module with flip-sign method
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalBasisFlipSign1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalTensorBasisFlipSign2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalTensorBasisFlipSign3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalTensorBasisFlipSign4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign5d"
   elseif (self._ndim == 6) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign6d"
   end

   self._flipSign = _m[self._polyOrder] -- function to flip sign   
end

function CartModalTensor:id() return "tensor" end
function CartModalTensor:ndim() return self._ndim end
function CartModalTensor:polyOrder() return self._polyOrder end
function CartModalTensor:numBasis() return self._numBasis end
function CartModalTensor:numSurfBasis() return self._numSurfBasis end
function CartModalTensor:evalBasis(z, b) return self._evalBasisFunc(z, b) end
function CartModalTensor:flipSign(dir, fIn, fOut) self._flipSign(dir, fIn, fOut) end

return {
   CartModalTensor = CartModalTensor   
}
