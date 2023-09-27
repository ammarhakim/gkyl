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
  GKYL_BC_GAUSSIAN = 1,
  GKYL_BC_MAXWELLIAN = 2};

enum gkyl_bc_emission_spectrum_gamma_type {
  GKYL_BC_FURMAN_PIVI = 0,
  GKYL_BC_SCHOU = 1,
  GKYL_BC_CONSTANT = 2};

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
  enum gkyl_bc_emission_spectrum_type bctype,
  enum gkyl_bc_emission_spectrum_gamma_type gammatype,
  double *bc_param, double *sey_param, int cdim, int vdim, bool use_gpu);

void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *f_skin, const struct gkyl_array *f_proj, struct gkyl_array *f_buff,
  struct gkyl_array *weight, struct gkyl_array *k,
  const struct gkyl_array *flux, struct gkyl_rect_grid *grid, struct gkyl_array *gamma,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_range *conf_r,
  const struct gkyl_range *buff_r);

void gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *gamma, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r);

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
   self._cDim = assert(tbl.cdim,
      "Updater.EmissionSpectrumBc: Must specify configuration space dimensions with 'cdim'.")
   local vDim = assert(tbl.vdim,
      "Updater.EmissionSpectrumBc: Must specify velocity space dimensions with 'vdim'.")
   assert(self._edge == "lower" or self._edge == "upper",
      "Updater.EmissionSpectrumBc: 'edge' must be 'lower' or 'upper'.")

   self._bcType  = assert(tbl.bcType,
      "Updater.EmissionSpectrumBc: Must specify BC type in 'bcType'.")
   self._gammaType  = assert(tbl.gammaType,
			     "Updater.EmissionSpectrumBc: Must specify SEY model type in 'gammaType'.")
   self._bcParam  = assert(tbl.bcParam,
      "Updater.EmissionSpectrumBc: Must specify BC type in 'bcType'.")
   self._gammaParam  = assert(tbl.gammaParam,
			     "Updater.EmissionSpectrumBc: Must specify SEY model type in 'gammaType'.")
   local onGrid = assert(tbl.onGrid,
      "Updater.EmissionSpectrumBc: Must specify the field we'll apply BCs to in 'onGrid'.")
   local onField = assert(tbl.onField,
      "Updater.EmissionSpectrumBc: Must specify the field we'll apply BCs to in 'onField'.")

   local edge = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local otherRange = onField:localRange()

   local bctype -- Match gkyl_bc_emission_spectrum_type in gkylzero/zero/gkyl_bc_emission_spectrum.h
       if self._bcType == "chung-everhart" then bctype = 0
   elseif self._bcType == "gaussian" then bctype = 1
   elseif self._bcType == "maxwellian" then bctype = 2
   end
   local gammatype -- Match gkyl_bc_emission_spectrum_gamma_type in gkylzero/zero/gkyl_bc_emission_spectrum.h
       if self._gammaType == "furman-pivi" then gammatype = 0
   elseif self._gammaType == "schou" then gammatype = 1
   elseif self._gammaType == "constant" then gammatype = 2
   end

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)
   
   self._zero = ffi.gc(
      ffiC.gkyl_bc_emission_spectrum_new(self._dir-1, edge, bctype, gammatype, self._bcParam, self._gammaParam, self._cDim, vDim, useGPU or 0),
      ffiC.gkyl_bc_emission_spectrum_release
   )
   if useGPU then
      ffiC.gkyl_bc_emission_spectrum_sey_calc(self._zero, onField._zeroDevice, onGrid._zero, otherRange)
   else
      ffiC.gkyl_bc_emission_spectrum_sey_calc(self._zero, onField._zero, onGrid._zero, otherRange)
   end
   
   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function EmissionSpectrumBc:_advance(tCurr, inFld, outFld)
   local fOther = assert(inFld[1],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local param = assert(inFld[2],
      "EmissionSpectrumBc.advance: Must-specify BC params.")
   local fProj = assert(inFld[3],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local flux = assert(inFld[4],
      "EmissionSpectrumBc.advance: Must-specify flux.")
   local inGrid = assert(inFld[5],
      "EmissionSpectrumBc.advance: Must-specify grid.")
   local gamma = assert(inFld[6],
      "EmissionSpectrumBc.advance: Must-specify gamma.")
   local weight     = assert(outFld[1],
      "EmissionSpectrumBc.advance: Must-specify weight output field")
   local k     = assert(outFld[2],
      "EmissionSpectrumBc.advance: Must-specify k output field")
   local fBuff     = assert(outFld[3],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   local otherRange = assert(inFld[7],
      "EmissionSpectrumBc.advance: Must-specify other species range")
   local boundRange = gamma:localRange()
   local confRange = weight:localRange()
   local buffRange = fBuff:localRange()
   ffiC.gkyl_bc_emission_spectrum_advance(self._zero, fOther._zero, fProj._zero, fBuff._zero, weight._zero, k._zero, flux._zero, inGrid._zero, gamma._zero, otherRange, boundRange, confRange, buffRange)
end

function EmissionSpectrumBc:_advanceOnDevice(tCurr, inFld, outFld)
   local fOther = assert(inFld[1],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local param = assert(inFld[2],
      "EmissionSpectrumBc.advance: Must-specify BC params.")
   local fProj = assert(inFld[3],
      "EmissionSpectrumBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local flux = assert(inFld[4],
      "EmissionSpectrumBc.advance: Must-specify flux.")
   local inGrid = assert(inFld[5],
      "EmissionSpectrumBc.advance: Must-specify grid.")
   local gamma = assert(inFld[6],
      "EmissionSpectrumBc.advance: Must-specify gamma.")
   local weight     = assert(outFld[1],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   local k     = assert(outFld[2],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   local fBuff     = assert(outFld[3],
      "EmissionSpectrumBc.advance: Must-specify an output field")
   local otherRange = assert(inFld[7],
      "EmissionSpectrumBc.advance: Must-specify other species range")
   local boundRange = gamma:localRange()
   local confRange = weight:localRange()
   local buffRange = fBuff:localRange()
   ffiC.gkyl_bc_emission_spectrum_advance(self._zero, fOther._zeroDevice, fProj._zeroDevice, fBuff._zeroDevice, weight._zeroDevice, k._zeroDevice, flux._zeroDevice, inGrid._zero, gamma._zeroDevice, otherRange, boundRange, confRange, buffRange)
end

function EmissionSpectrumBc:getDir() return self._dir end

function EmissionSpectrumBc:getEdge() return self._edge end

function EmissionSpectrumBc:label() return "Flux"..self._dirlabel..self._edge end

return EmissionSpectrumBc
