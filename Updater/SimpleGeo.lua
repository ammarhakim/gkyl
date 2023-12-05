-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local DataStruct  = require "DataStruct"
local ffi         = require "ffi"
local Lin          = require "Lib.Linalg"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
typedef struct gk_geometry gk_geometry;

/**
 * Create a new wave geometry object. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param range_ext 
 * @param basis configuration space basis
 * @param mapc2p Mapping from computational to physical space
 * @param mapc2p_ctx Context for use in mapping
 * @param bmag function which gives |B| in computational space
 * @param bmag_ctx Context for calculating |B|
 */
struct gk_geometry* gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool use_gpu);

/**
 * Create a new wave geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_wave_geom_release method.
 *
 * @param up Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gk_geometry* gkyl_gk_geometry_acquire(const struct gk_geometry* up);


void gkyl_gk_geometry_free(const struct gkyl_ref_count *ref);

/**
 * Release geometry object.
 *
 * @param wg Wave geometry object to release.
 */
void gkyl_gk_geometry_release(const struct gk_geometry *up);



void gkyl_gk_geometry_advance(struct gk_geometry* up, 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p, 
  struct gkyl_array *bmag, struct gkyl_array *g_ij, 
  struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_inv, 
  struct gkyl_array *gij, struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot, 
  struct gkyl_array *jacobtot_inv, struct gkyl_array *bmag_inv, struct gkyl_array *bmag_inv_sq, 
  struct gkyl_array *gxxj, struct gkyl_array *gxyj, struct gkyl_array *gyyj) ;

const struct gkyl_range* gkyl_gk_geometry_get_nrange(struct gk_geometry* geo){
  return geo->nrange;
}

struct gkyl_array* gkyl_gk_geometry_get_mc2p_nodal_fd(struct gk_geometry* geo){
  return geo->mc2p_nodal_fd;
}

double* gkyl_gk_geometry_get_dzc(struct gk_geometry* geo){
  return geo->dzc;
}
]]

-- Boundary condition updater.
local SimpleGeo = Proto(UpdaterBase)

function SimpleGeo:init(tbl)
   SimpleGeo.super.init(self, tbl) -- Setup base object.
   self.grid = tbl.grid
   self.basis = tbl.basis
   self.localRange = tbl.localRange
   self.localRangeExt = tbl.localRangeExt

   -- Set up BCs 0=periodic, 1=nonperiodic
   -- BCs are used to set ranges for calculating geometry
   -- and to choose stencils when calculating metrics
   -- Currently we have not thought of a case where we need to calculate geometry with periodic bcs
   self.BCs = Lin.IntVec(self.grid:ndim())
   for i = 1, self.grid:ndim() do
      self.BCs[i] = 1
   end
   self.BCs[2] = 0
   self.conversionRange = self.localRange
   
   print("ok, calling zer new for simple geo")
   self._zero = ffi.gc(ffiC.gkyl_gk_geometry_new(self.grid._zero, self.localRange, self.localRangeExt, self.basis._zero, nil, nil, nil, nil, false), ffiC.gkyl_gk_geometry_release)

   print("Finished init method for SimpleGeo")



end


function SimpleGeo:_advance(tCurr, inFlds, outFlds)
   function fib(n)
      if n == 1 or n == 2 then
         return 1,1
      end
      prev, prevPrev = fib(n-1)
      return prev+prevPrev, prev
   end
   local mapc2p_field, gFld, jacobGeo, jacobGeoInv, jacobTot, jacobTotInv, bmagInv, bmagInvSq, gxxJ, gxyJ, gyyJ, grFld, b_i, cmag, bmag = outFlds[1], outFlds[2], outFlds[3], outFlds[4], outFlds[5], outFlds[6], outFlds[7], outFlds[8], outFlds[9], outFlds[10], outFlds[11], outFlds[12], outFlds[13], outFlds[14], outFlds[15]

   ffiC.gkyl_gk_geometry_advance( self._zero, nil, nil, nil, nil, mapc2p_field._zero, bmag._zero, gFld._zero, jacobGeo._zero, jacobGeoInv._zero, grFld._zero, b_i._zero, cmag._zero, jacobTot._zero, jacobTotInv._zero, bmagInv._zero, bmagInvSq._zero, gxxJ._zero, gxyJ._zero, gyyJ._zero);

end

return SimpleGeo
 
