-- Gkyl ------------------------------------------------------------------------
--
-- Apply a sheath boundary condition in the parallel direction in
-- a gyrokinetic model.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local ffi            = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_bc_sheath_gyrokinetic gkyl_bc_sheath_gyrokinetic;

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
struct gkyl_bc_sheath_gyrokinetic* gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim, double q2Dm, bool use_gpu);

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
void gkyl_bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up);
]]

-- Boundary condition updater.
local GkSheathBc = Proto(UpdaterBase)

function GkSheathBc:init(tbl)
   GkSheathBc.super.init(self, tbl) -- Setup base object.

   self._grid    = assert(tbl.onGrid, "Updater.GkSheathBc: Must specify grid to use with 'onGrid'.")
   self._dir     = assert(tbl.dir, "Updater.GkSheathBc: Must specify direction to apply BCs with 'dir'.")
   self._edge    = assert(tbl.edge, "Updater.GkSheathBc: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   self._basis   = assert(tbl.basis, "Updater.GkSheathBc: Must specify the basis in 'basis'.")
   self._phiWall = assert(tbl.phiWall, "Updater.GkSheathBc: Must specify the wall potential in 'phiWall'.")
   local cDim    = assert(tbl.cdim, "Updater.GkSheathBc: Must specify configuration space dimensions with 'cdim'.")
   local onField = assert(tbl.onField, "Updater.GkSheathBc: Must specify the field we'll apply BCs to in 'onField'.")
   local mass    = assert(tbl.mass, "Updater.GkSheathBc: Must specify the species mass with 'mass'.")
   local charge  = assert(tbl.charge, "Updater.GkSheathBc: Must specify the species charge with 'charge'.")
   local useGPU  = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or 0)

   assert(self._edge == "lower" or self._edge == "upper", "Updater.GkSheathBc: 'edge' must be 'lower' or 'upper'.")

   local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localExtRange = onField:localExtRange()
   local numGhostVec   = self._edge == 'lower' and onField:lowerGhostVec() or onField:upperGhostVec()

   local q2Dm  = 2.*charge/mass
   local basis = useGPU and self._basis._zeroDevice or self._basis._zero

   self._zero = ffi.gc(ffiC.gkyl_bc_sheath_gyrokinetic_new(self._dir-1, edge, localExtRange, numGhostVec:data(),
                                                           basis, self._grid._zero, cDim, q2Dm, useGPU),
                       ffiC.gkyl_bc_sheath_gyrokinetic_release)

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function GkSheathBc:_advance(tCurr, inFld, outFld)
   local phi  = assert(inFld[1], "GkSheathBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut = assert(outFld[1], "GkSheathBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_sheath_gyrokinetic_advance(self._zero, phi._zero, self._phiWall._zero, qOut._zero)
end

function GkSheathBc:_advanceOnDevice(tCurr, inFld, outFld)
   local phi  = assert(inFld[1], "GkSheathBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut = assert(outFld[1], "GkSheathBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_sheath_gyrokinetic_advance(self._zero, phi._zeroDevice, self._phiWall._zeroDevice, qOut._zeroDevice)
end

function GkSheathBc:getDir() return self._dir end

function GkSheathBc:getEdge() return self._edge end

function GkSheathBc:label() return "Flux"..self._dirlabel..self._edge end

return GkSheathBc