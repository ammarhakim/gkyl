-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local DataStruct  = require "DataStruct"
local ffi         = require "ffi"
local Lin          = require "Lib.Linalg"
local Grid = require "Grid"

local ffiC = ffi.C
require "Lib.ZeroUtil"

local new, sizeof, typeof, metatype = xsys.from(ffi, "new, sizeof, typeof, metatype")


-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_efit gkyl_efit;

struct gkyl_efit{
  const struct gkyl_basis *rzbasis;
  bool use_gpu;
  char* filepath;
  int nr, nz;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  double rmin, rmax, zmin, zmax;
  double *rzlower;
  double *rzupper;
  int *rzcells;
  int *rzghost;
};


/**
 * Create new updater to project psi, psi/R, psi/R^2
 * new method fills info to contruct a grid
 *
 * @param rzbasis Basis object 
 * @param rz grid to be filled from efit
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_efit* gkyl_efit_new(const char *filepath, const struct gkyl_basis *rzbasis, bool use_gpu);



/**
 * Project psi, psi/R, psi/R^2
 *
 * @param psi gkyl_array to be filled
 */

void gkyl_efit_advance(gkyl_efit* up, struct gkyl_rect_grid* rzgrid, struct gkyl_range* rzlocal, struct gkyl_range* rzlocal_ext, struct gkyl_array* psizr, struct gkyl_array* psibyrzr,struct gkyl_array* psibyr2zr);

void gkyl_efit_release(gkyl_efit* up);

]]

-- Boundary condition updater.
local Efit = Proto(UpdaterBase)


local function createField(grid, basis, ghostCells, vComp, periodicSync, useDevice)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid           = grid,
      numComponents    = basis:numBasis()*vComp,
      ghost            = ghostCells,
      metaData         = {polyOrder = basis:polyOrder(),
                          basisType = basis:id()},
      syncPeriodicDirs = periodicSync,
      useDevice        = useDevice,
   }
   fld:clear(0.0)
   return fld
end

function Efit:init(tbl)
   Efit.super.init(self, tbl) -- Setup base object.
   self.efitFile = tbl.efitFile
   self.rzBasis = tbl.rzBasis

   self._zero = ffi.gc(ffiC.gkyl_efit_new(self.efitFile, self.rzBasis._zero, false), ffiC.gkyl_efit_release)

   self.rzCells = {}
   self.rzLower = {}
   self.rzUpper = {}
   self.rzGhost = {}

   for i = 1, 2 do
      self.rzCells[i] = self._zero.rzcells[i-1]
      self.rzLower[i] = self._zero.rzlower[i-1]
      self.rzUpper[i] = self._zero.rzupper[i-1]
   end

   for i = 1, 2 do
      self.rzGhost[i] = self._zero.rzghost[i-1]
   end

   self.rzGrid = Grid.RectCart{
      lower = self.rzLower,
      upper = self.rzUpper,
      cells = self.rzCells,
   }

   local lower = self.rzGrid._lower
   local upper = self.rzGrid._upper


   local psiRZ = createField(self.rzGrid, self.rzBasis, self.rzGhost, 1, false)
   self.rzLocalRange = psiRZ:localRange()
   self.rzLocalRangeExt = psiRZ:localExtRange()

   self.B0 = self._zero.bcentr
   self.R0 = self._zero.rcentr

end

function Efit:getRZGrid()
   return self.rzGrid
end

function Efit:bphifunc(t,xn)
   local R = xn[1]
   return self.B0*self.R0/R
end

function Efit:_advance(tCurr, inFlds, outFlds)
   local psiRZ, psibyrRZ, psibyr2RZ = outFlds[1], outFlds[2], outFlds[3]
   ffiC.gkyl_efit_advance( self._zero, self.rzGrid._zero, self.rzLocalRange, self.rzLocalRangeExt, psiRZ._zero, psibyrRZ._zero, psibyr2RZ._zero)
end

return Efit
 
