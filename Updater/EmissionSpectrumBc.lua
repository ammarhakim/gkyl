-- Gkyl ------------------------------------------------------------------------
--
-- Apply a basic boundary condition, in which the function in the ghost cell
-- is only a function of the skin cell next to it.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Grid        = require "Grid"
local Lin         = require "Lib.Linalg"
local Proto       = require "Lib.Proto"
local Range       = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// BC types in this updater.
enum gkyl_bc_emission_spectrum_type {
  GKYL_BC_CHUNG_EVERHART = 0,
  GKYL_BC_GAUSSIAN = 1};

// Object type
typedef struct gkyl_bc_emission_spectrum gkyl_bc_emission_spectrum;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param grid cartesian grid dynamic field is defined on.
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc).
 * @param local_conf_range_ext Local extended configuration range.
 * @param local_range_ext Local extended phase range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param bctype BC type (see gkyl_bc_emission_spectrum_type).
 * @param cbasis Configuration basis on which coefficients in array are expanded
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param cdim Configuration space dimensions.
 * @param vdim Velocity space dimensions.
 * @param bc_param Parameters used for calculating BCs.
 * @param gain Array of secondary electron gain values at the cell centers.
 * @param elastic Array of elastic backscattering gain values at the cell centers.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_emission_spectrum* gkyl_bc_emission_spectrum_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *ghost_r,
  enum gkyl_bc_emission_spectrum_type bctype, const struct gkyl_basis *cbasis, const struct gkyl_basis *basis,
  int cdim, int vdim, bool use_gpu);

double gkyl_bc_emission_spectrum_advance_cross(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *f_self, struct gkyl_array *f_other, struct gkyl_array *f_proj,
  double *bc_param, double flux, struct gkyl_rect_grid *grid, double *gain, const struct gkyl_range *other_r);

/**
 * Advance boundary conditions. Fill buffer array based on boundary conditions and copy
 * contents to ghost cells of input f_arr
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr, struct gkyl_array *f_proj, double k);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
]]

-- Boundary condition updater.
local EmissionSpectrumBc = Proto(UpdaterBase)

function EmissionSpectrumBc:init(tbl)
   EmissionSpectrumBc.super.init(self, tbl) -- Setup base object.
   self._dir  = assert(tbl.dir,
      "Updater.EmissionSpectrumBc: Must specify direction to apply BCs with 'dir'.")
   self._edge = assert(tbl.edge,
      "Updater.EmissionSpectrumBc: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   local cDim = assert(tbl.cdim,
      "Updater.EmissionSpectrumBc: Must specify configuration space dimensions with 'cdim'.")
   local vDim = assert(tbl.vdim,
      "Updater.EmissionSpectrumBc: Must specify velocity space dimensions with 'vdim'.")
   assert(self._edge == "lower" or self._edge == "upper",
      "Updater.EmissionSpectrumBc: 'edge' must be 'lower' or 'upper'.")

   self._bcType  = assert(tbl.bcType,
      "Updater.EmissionSpectrumBc: Must specify BC type in 'bcType'.")
   self._basis   = assert(tbl.basis,
      "Updater.EmissionSpectrumBc: Must specify the basis in 'basis'.")
   self._confBasis   = assert(tbl.confBasis,
      "Updater.EmissionSpectrumBc: Must specify the conf basis in 'confBasis'.")
   local onField = assert(tbl.onField,
      "Updater.EmissionSpectrumBc: Must specify the field we'll apply BCs to in 'onField'.")

   local edge = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localGhostRange = self._edge == 'lower' and onField:localGhostRangeLower()[self._dir] or onField:localGhostRangeUpper()[self._dir]

   local bctype -- Match gkyl_bc_emission_type in gkylzero/zero/gkyl_bc_emission.h
       if self._bcType == "chung" then bctype = 0
   elseif self._bcType == "gauss" then bctype = 1
   end

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)
   local basis  = useGPU and self._basis._zeroDevice or self._basis._zero
   local cbasis  = useGPU and self._confBasis._zeroDevice or self._confBasis._zero
   
   self._zero = ffi.gc(
      ffiC.gkyl_bc_emission_spectrum_new(self._dir-1, edge,
	localGhostRange, bctype, cbasis, basis,
	cDim, vDim, useGPU or 0),
      ffiC.gkyl_bc_emission_spectrum_release
   )

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function EmissionSpectrumBc:_advance(tCurr, inFld, outFld)
   local fOther = assert(inFld[1],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local fProj = assert(inFld[2],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local param = assert(inFld[3],
      "EmissionSpectrumBc.advance: Must-specify BC params.")
   local flux = assert(inFld[4],
      "EmissionSpectrumBc.advance: Must-specify flux.")
   local inGrid = assert(inFld[5],
      "EmissionSpectrumBc.advance: Must-specify grid.")
   local gain = assert(inFld[6],
      "EmissionSpectrumBc.advance: Must-specify gamma.")
   local qOut     = assert(outFld[1],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   local otherRange = self._edge == 'lower' and fOther:localSkinRangeLower()[self._dir] or fOther:localSkinRangeUpper()[self._dir]
   local k = ffiC.gkyl_bc_emission_spectrum_advance_cross(self._zero, qOut._zero, fOther._zero, fProj._zero, param:data(), flux, inGrid._zero, gain, otherRange)
   return k
end

-- function EmissionSpectrumBc:_advance(tCurr, inFld, outFld)
--    local bufferIn = assert(inFld[1],
--       "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
--    local qOut     = assert(outFld[1],
--       "EmissionSpectrumBc.advance: Must-specify an output field")
--    ffiC.gkyl_bc_emission_spectrum_advance(self._zero, bufferIn._zero, qOut._zero)
-- end

function EmissionSpectrumBc:_advanceOnDevice(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut = assert(outFld[1],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_emission_spectrum_advance(self._zero, bufferIn._zeroDevice, qOut._zeroDevice)
end

function EmissionSpectrumBc:getDir() return self._dir end

function EmissionSpectrumBc:getEdge() return self._edge end

function EmissionSpectrumBc:label() return "Flux"..self._dirlabel..self._edge end

return EmissionSpectrumBc
