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
enum gkyl_bc_see_type { GKYL_BC_SEE = 0 };

// Object type
typedef struct gkyl_bc_see gkyl_bc_see;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param grid cartesian grid dynamic field is defined on.
 * @param cdim Configuration space dimensions.
 * @param q2Dm charge-to-mass ratio times 2.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_see* gkyl_bc_see_new(const struct gkyl_rect_grid *grid, int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_range_ext, const struct gkyl_range *local_range_ext,
		const int *num_ghosts, const struct gkyl_basis *cbasis, const struct gkyl_basis *basis,
		int num_comp, int cdim, int vdim, const double *bc_param, double *gain, double *elastic, bool use_gpu);

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 */
void gkyl_bc_see_advance(const struct gkyl_bc_see *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_see_release(struct gkyl_bc_see *up);
]]

-- Boundary condition updater.
local SEEBc = Proto(UpdaterBase)

function SEEBc:init(tbl)
   SEEBc.super.init(self, tbl) -- Setup base object.

   self._grid = assert(tbl.onGrid, "Updater.SEEBc: Must specify grid to use with 'onGrid'.")
   self._dir  = assert(tbl.dir, "Updater.SEEBc: Must specify direction to apply BCs with 'dir'.")
   self._edge = assert(tbl.edge, "Updater.SEEBc: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   local cDim = assert(tbl.cdim, "Updater.SEEBc: Must specify configuration space dimensions with 'cdim'.")
   local vDim = assert(tbl.vdim, "Updater.SEEBc: Must specify configuration space dimensions with 'cdim'.")
   local gain = assert(tbl.gain, "Updater.SEEBc: Must specify configuration space dimensions with 'cdim'.")
   local elastic = tbl.elastic or nil
   assert(self._edge == "lower" or self._edge == "upper", "Updater.SEEBc: 'edge' must be 'lower' or 'upper'.")

   self._bcType  = assert(tbl.bcType, "Updater.SEEBc: Must specify BC type in 'bcType'.")
   self._basis   = assert(tbl.basis, "Updater.SEEBc: Must specify the basis in 'basis'.")
   self._confBasis   = assert(tbl.confBasis, "Updater.SEEBc: Must specify the conf basis in 'confBasis'.")
   local onField = assert(tbl.onField, "Updater.SEEBc: Must specify the field we'll apply BCs to in 'onField'.")

   local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localExtRange = onField:localExtRange()
   local confLocalExtRange = localExtRange:selectFirst(cDim)
   local numGhostVec   = self._edge == 'lower' and onField:lowerGhostVec() or onField:upperGhostVec()

   local bctype -- Match gkyl_bc_emission_type in gkylzero/zero/gkyl_bc_emission.h
       if self._bcType == "gain"        then bctype = 0
   end

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)
   local basis  = useGPU and self._basis._zeroDevice or self._basis._zero
   local cbasis  = useGPU and self._confBasis._zeroDevice or self._confBasis._zero
   local bcParam = tbl.bcParam
   
   self._zero = ffi.gc(ffiC.gkyl_bc_see_new(self._grid._zero, self._dir-1, edge, confLocalExtRange, localExtRange, numGhostVec:data(), cbasis, basis, onField:numComponents(), cDim, vDim, bcParam:data(), gain, elastic, useGPU or 0),
                       ffiC.gkyl_bc_see_release)

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function SEEBc:_advance(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "SEEBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "SEEBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_see_advance(self._zero, bufferIn._zero, qOut._zero)
end

function SEEBc:_advanceOnDevice(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "SEEBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "SEEBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_see_advance(self._zero, bufferIn._zeroDevice, qOut._zeroDevice)
end

function SEEBc:getDir() return self._dir end

function SEEBc:getEdge() return self._edge end

function SEEBc:label() return "Flux"..self._dirlabel..self._edge end

return SEEBc
