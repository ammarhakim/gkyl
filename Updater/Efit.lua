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
  const char* filepath;
  int nr, nz;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  double rmin, rmax, zmin, zmax;
  double *rzlower;
  double *rzupper;
  int *rzcells;
  int *rzghost;

  const struct gkyl_basis *fluxbasis;
  double *fluxlower;
  double *fluxupper;
  int *fluxcells;
  int *fluxghost;
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
gkyl_efit* gkyl_efit_new(const char *filepath, const struct gkyl_basis *rzbasis, const struct gkyl_basis *fluxbasis, bool use_gpu);



/**
 * Project psi, psi/R, psi/R^2
 *
 * @param rzbasis Basis object 
 * @param rz grid to be filled from efit
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */

void gkyl_efit_advance(gkyl_efit* up, struct gkyl_rect_grid* rzgrid, struct gkyl_rect_grid* fluxgrid, struct gkyl_range* rzlocal, struct gkyl_range* rzlocal_ext, struct gkyl_array* psizr, struct gkyl_array* psibyrzr,struct gkyl_array* psibyr2zr, struct gkyl_range* fluxlocal, struct gkyl_range* fluxlocal_ext, struct gkyl_array* fpolflux, struct gkyl_array* qflux);

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
   self.fBasis = tbl.fBasis

   self._zero = ffi.gc(ffiC.gkyl_efit_new(self.efitFile, self.rzBasis._zero, self.fBasis._zero, false), ffiC.gkyl_efit_release)

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

   local psiRZ = createField(self.rzGrid, self.rzBasis, self.rzGhost, 1, false)
   self.rzLocalRange = psiRZ:localRange()
   self.rzLocalRangeExt = psiRZ:localExtRange()

   self.fCells = {}
   self.fLower = {}
   self.fUpper = {}
   self.fGhost = {}

   self.fCells[1] = self._zero.fluxcells[0]
   self.fLower[1] = self._zero.fluxlower[0]
   self.fUpper[1] = self._zero.fluxupper[0]
   self.fGhost[1] = self._zero.fluxghost[0]
   self.fGhost[2] = self._zero.fluxghost[1]

   self.fGrid = Grid.RectCart{
      lower = self.fLower,
      upper = self.fUpper,
      cells = self.fCells,
   }

   local fpol = createField(self.fGrid, self.fBasis, self.fGhost, 1, false)
   self.fLocalRange = fpol:localRange()
   self.fLocalRangeExt = fpol:localExtRange()

   self.B0 = self._zero.bcentr
   self.R0 = self._zero.rcentr
   self.psiSep = self._zero.sibry

end

function Efit:getRZGrid()
   return self.rzGrid
end

function Efit:getfGrid()
   return self.fGrid
end

function Efit:bphifunc(t,xn)
   local R = xn[1]
   return self.B0*self.R0/R
end

function Efit:_advance(tCurr, inFlds, outFlds)
   local psiRZ, psibyrRZ, psibyr2RZ, fpol, qflux = outFlds[1], outFlds[2], outFlds[3], outFlds[4], outFlds[5]
   ffiC.gkyl_efit_advance( self._zero, self.rzGrid._zero, self.fGrid._zero, self.rzLocalRange, self.rzLocalRangeExt, psiRZ._zero, psibyrRZ._zero, psibyr2RZ._zero, self.fLocalRange, self.fLocalRangeExt, fpol._zero, qflux._zero)
end

return Efit
 
