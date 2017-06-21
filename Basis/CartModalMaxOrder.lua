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

local CartModalMaxOrder = {}
function CartModalMaxOrder:new(tbl)
   local self = setmetatable({}, CartModalMaxOrder)

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

   _m = nil -- to store module with vol -> surf projections
   if (self._ndim == 1) then
      _m = require "Basis._data.ModalMaxOrderBasisSurf1d"
   elseif (self._ndim == 2) then
      _m = require "Basis._data.ModalMaxOrderBasisSurf2d"
   elseif (self._ndim == 3) then
      _m = require "Basis._data.ModalMaxOrderBasisSurf3d"
   elseif (self._ndim == 4) then
      _m = require "Basis._data.ModalMaxOrderBasisSurf4d"
   elseif (self._ndim == 5) then
      _m = require "Basis._data.ModalMaxOrderBasisSurf5d"
   elseif (self._ndim == 6) then
      --assert(false, "NDIM 6 NYI!!")
   end

   self._projectVolToSurfLower = {} -- functions to project volume expansion of lower surface
   self._projectVolToSurfUpper = {} -- functions to project volume expansion of upper surface
   for d = 1, self._ndim do
      self._projectVolToSurfLower[d] = _m[self._polyOrder][d].lower
      self._projectVolToSurfUpper[d] = _m[self._polyOrder][d].upper
   end   
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartModalMaxOrder, { __call = function (self, o) return self.new(self, o) end })

CartModalMaxOrder.__index = {
   id = function (self)
      return "maximal-order"
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
   volumeToLowerSurfExpansion = function (self, dir, volIn, surfOut)
      return self._projectVolToSurfLower[dir](volIn, surfOut)
   end,
   volumeToUpperSurfExpansion = function (self, dir, volIn, surfOut)
      return self._projectVolToSurfUpper[dir](volIn, surfOut)
   end,
}

return {
   CartModalMaxOrder = CartModalMaxOrder   
}
