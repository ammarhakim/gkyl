-- Gkyl ------------------------------------------------------------------------
--
-- Modal Serendipity elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"
local ffi   = require "ffi"
local ffiC  = ffi.C
require "Basis.BasisCdef"

ffi.cdef [[ 
/**
 * Create new modal serendipity basis function object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim, int poly_order);

/**
 * Create new hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_gk_hybrid(struct gkyl_basis *basis, int ndim);
]]

-- CartModalSerendipity -----------------------------------------------------------
--
-- Modal serendipity basis set. See: Found. Comput. Math 2011 11:337 Arnold & Awanou
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local CartModalSerendipity = Proto()
function CartModalSerendipity:init(tbl)
   -- read data from input table
   self._ndim = assert(tbl.ndim, "Basis.CartModalSerendipity: Must specify dimension using 'ndim'")
   self._polyOrder = assert(tbl.polyOrder, "Basis.CartModalSerendipity: Must specify polynonial order with 'polyOrder'")

   if (self._polyOrder < 0) or (self._polyOrder > 3) then
      assert(false, "Basis.CartModalSerendipity: Polynomial order must be between 0 and 3")
   end

   -- Option to use hybrid basis. It has vpar^2 & its tensor product w/ other monomials.
   self._isGkHybrid = tbl.useGkHybrid 

   -- create gkylzero gkyl_basis struct
   self._zero = ffi.new("struct gkyl_basis")

   if self._isGkHybrid then
      assert(self._polyOrder == 1, "Basis.CartModalSerendipity: GK hybrid basis requires polyOrder=1.")
      ffiC.gkyl_cart_modal_gk_hybrid(self._zero, self._ndim)
   else
      ffiC.gkyl_cart_modal_serendip(self._zero, self._ndim, self._polyOrder)
   end

   self._numBasis = self._zero.num_basis
end

function CartModalSerendipity:id() return self._isGkHybrid and "gk_hybrid" or "serendipity" end
function CartModalSerendipity:ndim() return self._ndim end
function CartModalSerendipity:polyOrder() return self._polyOrder end
function CartModalSerendipity:numBasis() return self._numBasis end
-- These three functions expect C pointers and a 1-indexed direction as input.
function CartModalSerendipity:evalBasis(z, b) return self._zero.eval(z, b) end
function CartModalSerendipity:flipSign(dir, fIn, fOut) self._zero.flip_odd_sign(dir-1, fIn, fOut) end
function CartModalSerendipity:flipEvenSign(dir, fIn, fOut) self._zero.flip_even_sign(dir-1, fIn, fOut) end

return {
   CartModalSerendipity = CartModalSerendipity   
}
