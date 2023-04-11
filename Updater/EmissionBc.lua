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
enum gkyl_bc_emission_type { GKYL_BC_CONSTANT_GAIN = 0 };

// Object type
typedef struct gkyl_bc_emission gkyl_bc_emission;

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
struct gkyl_bc_emission* gkyl_bc_emission_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, enum gkyl_bc_emission_type bctype, const struct gkyl_basis *basis,
  int num_comp, int cdim, double *bc_param, const struct gkyl_array *bc_field, bool use_gpu);

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
void gkyl_bc_emission_advance(const struct gkyl_bc_emission *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_release(struct gkyl_bc_emission *up);
]]

-- Boundary condition updater.
local EmissionBc = Proto(UpdaterBase)

function EmissionBc:init(tbl)
   EmissionBc.super.init(self, tbl) -- Setup base object.

   self._grid = assert(tbl.onGrid, "Updater.EmissionBc: Must specify grid to use with 'onGrid'.")
   self._dir  = assert(tbl.dir, "Updater.EmissionBc: Must specify direction to apply BCs with 'dir'.")
   self._edge = assert(tbl.edge, "Updater.EmissionBc: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   local cDim = assert(tbl.cdim, "Updater.EmissionBc: Must specify configuration space dimensions with 'cdim'.")
   assert(self._edge == "lower" or self._edge == "upper", "Updater.EmissionBc: 'edge' must be 'lower' or 'upper'.")

   self._bcType  = assert(tbl.bcType, "Updater.EmissionBc: Must specify BC type in 'bcType'.")
   self._basis   = assert(tbl.basis, "Updater.EmissionBc: Must specify the basis in 'basis'.")
   local onField = assert(tbl.onField, "Updater.EmissionBc: Must specify the field we'll apply BCs to in 'onField'.")
   self._bcField = assert(tbl.onField, "Updater.EmissionBc: Must specify the BC emission parameters in 'bcField'.")

   local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localExtRange = onField:localExtRange()
   local numGhostVec   = self._edge == 'lower' and onField:lowerGhostVec() or onField:upperGhostVec()

   local bctype -- Match gkyl_bc_emission_type in gkylzero/zero/gkyl_bc_emission.h
       if self._bcType == "gain"        then bctype = 0
   end

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)
   local basis  = useGPU and self._basis._zeroDevice or self._basis._zero
   local bcParam = tbl.bcParam
   
   self._zero = ffi.gc(ffiC.gkyl_bc_emission_new(self._dir-1, edge, localExtRange, numGhostVec:data(), bctype,
						 basis, onField:numComponents(), cDim, bcParam:data(), nil, useGPU or 0),
                       ffiC.gkyl_bc_emission_release)

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function EmissionBc:_advance(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "EmissionBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "EmissionBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_emission_advance(self._zero, bufferIn._zero, qOut._zero)
end

function EmissionBc:_advanceOnDevice(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "EmissionBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "EmissionBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_emission_advance(self._zero, bufferIn._zeroDevice, qOut._zeroDevice)
end

function EmissionBc:getDir() return self._dir end

function EmissionBc:getEdge() return self._edge end

function EmissionBc:label() return "Flux"..self._dirlabel..self._edge end

return EmissionBc
