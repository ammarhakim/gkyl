-- Gkyl ------------------------------------------------------------------------
--
-- Modal Serendipity elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"

-- CartModalSerendipity -----------------------------------------------------------
--
-- Modal serendipity basis set. See: Found. Comput. Math 2011 11:337 Arnold & Awanou
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local function numBasis(ndim, polyOrder)
   local numBasis_d2 = {4, 8, 12, 17, 23, 30}
   local numBasis_d3 = {8, 20, 32, 50, 74}
   local numBasis_d4 = {16, 48, 80, 136}
   local numBasis_d5 = {32, 112, 192}
   local numBasis_d6 = {64, 256}

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
      elseif (ndim == 6) then
	 nbasis = numBasis_d6[polyOrder]
      end
   end
   return nbasis
end

local CartModalSerendipity = Proto()
function CartModalSerendipity:init(tbl)
   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalSerendipity: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalSerendipity: Must specify polynonial order with 'polyOrder'")

   self._numBasis = numBasis(self._ndim, self._polyOrder)
   self._numSurfBasis = numBasis(self._ndim-1, self._polyOrder)

   local _m = nil -- to store module with evaluation code
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      assert(self._polyOrder <= 7, "For 1D polynomial order must be either 1, 2, 3, 4, 5, 6, or 7")
      _m = require "Basis._data.ModalBasis1d"
   elseif (self._ndim == 2) then
      assert(self._polyOrder <= 6, "For 2D polynomial order must be either 1, 2, 3, 4, 5, or 6")
      _m = require "Basis._data.ModalSerendipBasis2d"
   elseif (self._ndim == 3) then
      assert(self._polyOrder <= 5, "For 3D polynomial order must be either 1, 2, 3, 4 or 5")
      _m = require "Basis._data.ModalSerendipBasis3d"
   elseif (self._ndim == 4) then
      assert(self._polyOrder <= 4, "For 4D polynomial order must be either 1, 2, 3, or 4")
      _m = require "Basis._data.ModalSerendipBasis4d"
   elseif (self._ndim == 5) then
      assert(self._polyOrder <= 3, "For 5D polynomial order must be either 1, 2, or 3")
      _m = require "Basis._data.ModalSerendipBasis5d"
   elseif (self._ndim == 6) then
      assert(self._polyOrder <= 2, "For 6D polynomial order must be either 1 or 2")
      _m = require "Basis._data.ModalSerendipBasis6d"

   end

   self._evalBasisFunc = _m[self._polyOrder] -- function to evaluate basis functions

   _m = nil -- to store module with flip-sign method
   -- get handle to function to compute basis functions at specified coordinates   
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalBasisFlipSign1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign5d"
   elseif (self._ndim == 6) then
      _m = require "Basis._data.ModalSerendipBasisFlipSign6d"
   end

   self._flipSign = function(dir, f, out) out[1] = f[1] end
   if self._numBasis > 1 then
      self._flipSign = _m[self._polyOrder] -- function to flip sign
   end
end

function CartModalSerendipity:id() return "serendipity" end
function CartModalSerendipity:ndim() return self._ndim end
function CartModalSerendipity:polyOrder() return self._polyOrder end
function CartModalSerendipity:numBasis() return self._numBasis end
function CartModalSerendipity:numSurfBasis() return self._numSurfBasis end
function CartModalSerendipity:evalBasis(z, b) return self._evalBasisFunc(z, b) end
function CartModalSerendipity:flipSign(dir, fIn, fOut) self._flipSign(dir, fIn, fOut) end

return {
   CartModalSerendipity = CartModalSerendipity   
}
