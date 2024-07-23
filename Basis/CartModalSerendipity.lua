-- Gkyl ------------------------------------------------------------------------
--
-- Modal Serendipity elements on Cartesian meshes.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local xsys = require "xsys"

-- Gkyl libraries
local Proto = require "Lib.Proto"
local ffi   = require "ffi"
local ffiC  = ffi.C
require "Basis.BasisCdef"

ffi.cdef [[ 
/**
 * Assign object members in modal serendipity basis object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 */
void gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim,
  int poly_order);
void gkyl_cart_modal_serendip_cu_dev(struct gkyl_basis *basis, int ndim,
  int poly_order);

/**
 * Create new modal serendipity basis function object.
 * This basis needs to be deallocated with free/release methods.
 *
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_serendip_new(int ndim, int poly_order);
struct gkyl_basis * gkyl_cart_modal_serendip_cu_dev_new(int ndim, int poly_order);

/**
 * Free the memory associated with basis objects
 * created with _new methods above.
 *
 * @param basis Basis object to free.
 */
void gkyl_cart_modal_basis_release(struct gkyl_basis *basis);
void gkyl_cart_modal_basis_release_cu(struct gkyl_basis *basis);
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

   -- create gkylzero gkyl_basis struct
   self._zero = ffi.new("struct gkyl_basis")
   ffiC.gkyl_cart_modal_serendip(self._zero, self._ndim, self._polyOrder)

   if xsys.pickBool(tbl.useDevice, GKYL_USE_GPU) then
      self._zeroDevice = ffi.gc(ffiC.gkyl_cart_modal_serendip_cu_dev_new(self._ndim, self._polyOrder),
                                ffiC.gkyl_cart_modal_basis_release_cu)
   else
      self._zeroDevice = nil
   end

   self._numBasis = self._zero.num_basis
end

function CartModalSerendipity:id() return "serendipity" end
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
