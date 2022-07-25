-- Gkyl ------------------------------------------------------------------------
--
-- Modal GkHybrid elements on Cartesian meshes. This is a basis consisting of
-- p=2 elements in vpar, and p=1 in other dimensions (x,y,z,mu).
--
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
 * Create new hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included.
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_gkhybrid(struct gkyl_basis *basis, int cdim, int vdim);
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
