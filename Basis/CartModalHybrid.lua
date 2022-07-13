-- Gkyl ------------------------------------------------------------------------
--
-- Modal Hybrid elements on Cartesian meshes. This is a basis consisting of
-- p=1 elements in configuration space, and p=2 Serendipity elements in
-- velocity space.
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
 * Create new hybrid basis. These are p=1 in configuration space
 * and p=2 in velocity space.
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_hybrid(struct gkyl_basis *basis, int cdim, int vdim);
]]

-- CartModalHybrid ----------------------------------------------------------------
--
-- Modal hybrid basis set. This is p=1 in position, p=2 Serendipity in velocity.
-----------------------------------------------------------------------------------

-- Number of basis function in ndim dimensions with specfied polyOrder
local CartModalHybrid = Proto()
function CartModalHybrid:init(tbl)
   -- read data from input table
   self._cdim = assert(tbl.cdim, "Basis.CartModalHybrid: Must specify configuration space dimension using 'cdim'")
   self._vdim = assert(tbl.vdim, "Basis.CartModalHybrid: Must specify velocity space dimension using 'vdim'")

   self._polyOrder = 1

   self._ndim = self._cdim + self._vdim

   -- create gkylzero gkyl_basis struct
   self._zero = ffi.new("struct gkyl_basis")

   ffiC.gkyl_cart_modal_hybrid(self._zero, self._cdim, self._vdim)

   self._numBasis = self._zero.num_basis
end

function CartModalHybrid:id() return "hybrid" end
function CartModalHybrid:ndim() return self._ndim end
function CartModalHybrid:polyOrder() return self._polyOrder end
function CartModalHybrid:numBasis() return self._numBasis end
-- These three functions expect C pointers and a 1-indexed direction as input.
function CartModalHybrid:evalBasis(z, b) return self._zero.eval(z, b) end
function CartModalHybrid:flipSign(dir, fIn, fOut) self._zero.flip_odd_sign(dir-1, fIn, fOut) end
function CartModalHybrid:flipEvenSign(dir, fIn, fOut) self._zero.flip_even_sign(dir-1, fIn, fOut) end

return {
   CartModalHybrid = CartModalHybrid   
}
