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

// Object type
typedef struct gkyl_deflate_geo gkyl_deflate_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */


gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis,const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, bool use_gpu);
/**
 * Advance deflate_geo (compute the derived_geo coefficients).
 *
 * @param up deflate_geo updater object.
 * @param crange Config-space range.
 * @param gFld field containing DG rep of the metric coefficients
 * @param jFld output field where jacobian will be placed
 */


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_release(gkyl_deflate_geo* up);

]]

-- Boundary condition updater.
local GeometryDeflate = Proto(UpdaterBase)

function GeometryDeflate:init(tbl)
   GeometryDeflate.super.init(self, tbl) -- Setup base object.

   self.grid = tbl.grid
   self.deflatedGrid = tbl.deflatedGrid
   self.basis = tbl.basis
   self.deflatedBasis = tbl.deflatedBasis
   self.localRange = tbl.localRange
   self.localRangeExt = tbl.localRangeExt
   self.deflatedLocalRange = tbl.deflatedLocalRange
   self.deflatedLocalRangeExt = tbl.deflatedLocalRangeExt
   self.remDirs = tbl.remDirs



   self._zero_deflate = ffi.gc(ffiC.gkyl_deflate_geo_new(self.basis._zero, self.deflatedBasis._zero, self.grid._zero, self.deflatedGrid._zero, self.remDirs:data(), false), ffiC.gkyl_deflate_geo_release)
   
   


end

function GeometryDeflate:_advance(tCurr, inFlds, outFlds)
   local ncomp = Lin.IntVec(19)
   for i = 1,#outFlds do
      ncomp[i] = inFlds[i]:numComponents()/self.basis:numBasis()
   end
   for i = 1, #outFlds do
      ffiC.gkyl_deflate_geo_advance(self._zero_deflate, self.localRange, self.deflatedLocalRange, inFlds[i]._zero, outFlds[i]._zero, ncomp:data()[i-1])
   end



end

return GeometryDeflate
 
