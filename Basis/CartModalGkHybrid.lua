-- Gkyl ------------------------------------------------------------------------
--
-- Modal GkHybrid elements on Cartesian meshes. This is a basis consisting of
-- p=2 elements in vpar, and p=1 in other dimensions (x,y,z,mu).
--
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
 * Assign object members in hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included.
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 */
void gkyl_cart_modal_gkhybrid(struct gkyl_basis *basis, int cdim, int vdim);
void gkyl_cart_modal_gkhybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim);

/**
 * Create new hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included.
 * This basis needs to be deallocated with free/release methods.
 *
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_gkhybrid_new(int cdim, int vdim);
struct gkyl_basis * gkyl_cart_modal_gkhybrid_cu_dev_new(int cdim, int vdim);

/**
 * Free the memory associated with basis objects
 * created with _new methods above.
 *
 * @param basis Basis object to free.
 */
void gkyl_cart_modal_basis_release(struct gkyl_basis *basis);
void gkyl_cart_modal_basis_release_cu(struct gkyl_basis *basis);
]]

-- CartModalGkHybrid --------------------------------------------------------------
--
-- Modal gkhybrid basis set. This is p=1 in vpar, p=2 in other dimensiona.
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local CartModalGkHybrid = Proto()
function CartModalGkHybrid:init(tbl)
   -- read data from input table
   self._cdim = assert(tbl.cdim, "Basis.CartModalGkHybrid: Must specify configuration space dimension using 'cdim'")
   self._vdim = assert(tbl.vdim, "Basis.CartModalGkHybrid: Must specify velocity space dimension using 'vdim'")

   self._polyOrder = 1

   self._ndim = self._cdim + self._vdim

   -- create gkylzero gkyl_basis struct
   self._zero = ffi.new("struct gkyl_basis")
   ffiC.gkyl_cart_modal_gkhybrid(self._zero, self._cdim, self._vdim)

   if xsys.pickBool(tbl.useDevice, GKYL_USE_GPU) then
      self._zeroDevice = ffi.gc(ffiC.gkyl_cart_modal_gkhybrid_cu_dev_new(self._cdim, self._vdim),
                                ffiC.gkyl_cart_modal_basis_release_cu)
   else
      self._zeroDevice = nil
   end

   self._numBasis = self._zero.num_basis
end

function CartModalGkHybrid:id() return "gkhybrid" end
function CartModalGkHybrid:ndim() return self._ndim end
function CartModalGkHybrid:polyOrder() return self._polyOrder end
function CartModalGkHybrid:numBasis() return self._numBasis end
-- These three functions expect C pointers and a 1-indexed direction as input.
function CartModalGkHybrid:evalBasis(z, b) return self._zero.eval(z, b) end
function CartModalGkHybrid:flipSign(dir, fIn, fOut) self._zero.flip_odd_sign(dir-1, fIn, fOut) end
function CartModalGkHybrid:flipEvenSign(dir, fIn, fOut) self._zero.flip_even_sign(dir-1, fIn, fOut) end

return {
   CartModalGkHybrid = CartModalGkHybrid   
}
