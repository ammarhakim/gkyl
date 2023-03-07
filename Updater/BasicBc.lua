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
// Constants to represent lower/upper edges
enum gkyl_edge_loc { GKYL_LOWER_EDGE = 0, GKYL_UPPER_EDGE = 1 };

// BC types in this updater.
enum gkyl_bc_basic_type { 
   GKYL_BC_COPY = 0, 
   GKYL_BC_ABSORB = 1, 
   GKYL_BC_REFLECT = 2,
   GKYL_BC_MAXWELL_PEC = 3, 
};

// Object type
typedef struct gkyl_bc_basic gkyl_bc_basic;

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param bctype BC type (see gkyl_bc_basic_type).
 * @param basis Basis on which coefficients in array are expanded.
 * @param num_comp Number of components (DOFs) within a cell.
 * @param cdim Configuration space dimensions.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_basic* gkyl_bc_basic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range* local_range_ext,
  const int *num_ghosts, enum gkyl_bc_basic_type bctype, const struct gkyl_basis *basis, int num_comp, int cdim, bool use_gpu);

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_basic_advance(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_basic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_basic_release(struct gkyl_bc_basic *up);
]]

-- Boundary condition updater.
local BasicBc = Proto(UpdaterBase)

function BasicBc:init(tbl)
   BasicBc.super.init(self, tbl) -- Setup base object.

   self._grid = assert(tbl.onGrid, "Updater.BasicBc: Must specify grid to use with 'onGrid'.")
   self._dir  = assert(tbl.dir, "Updater.BasicBc: Must specify direction to apply BCs with 'dir'.")
   self._edge = assert(tbl.edge, "Updater.BasicBc: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   local cDim = assert(tbl.cdim, "Updater.BasicBc: Must specify configuration space dimensions with 'cdim'.")
   assert(self._edge == "lower" or self._edge == "upper", "Updater.BasicBc: 'edge' must be 'lower' or 'upper'.")

   self._bcType  = assert(tbl.bcType, "Updater.BasicBc: Must specify BC type in 'bcType'.")
   self._basis   = assert(tbl.basis, "Updater.BasicBc: Must specify the basis in 'basis'.")
   local onField = assert(tbl.onField, "Updater.BasicBc: Must specify the field we'll apply BCs to in 'onField'.")

   local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localExtRange = onField:localExtRange()
   local numGhostVec   = self._edge == 'lower' and onField:lowerGhostVec() or onField:upperGhostVec()

   local bctype -- Match gkyl_bc_basic_type in gkylzero/zero/gkyl_bc_basic.h
       if self._bcType == "copy"        then bctype = 0
   elseif self._bcType == "absorb"      then bctype = 1 
   elseif self._bcType == "reflect"     then bctype = 2 
   elseif self._bcType == "maxwell_pec" then bctype = 3 
   end

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)
   local basis  = useGPU and self._basis._zeroDevice or self._basis._zero

   self._zero = ffi.gc(ffiC.gkyl_bc_basic_new(self._dir-1, edge, localExtRange, numGhostVec:data(), bctype,
                                              basis, onField:numComponents(), cDim, useGPU or 0),
                       ffiC.gkyl_bc_basic_release)

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function BasicBc:_advance(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "BasicBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "BasicBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_basic_advance(self._zero, bufferIn._zero, qOut._zero)
end

function BasicBc:_advanceOnDevice(tCurr, inFld, outFld)
   local bufferIn = assert(inFld[1], "BasicBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
   local qOut     = assert(outFld[1], "BasicBc.advance: Must-specify an output field")
   ffiC.gkyl_bc_basic_advance(self._zero, bufferIn._zeroDevice, qOut._zeroDevice)
end

function BasicBc:getDir() return self._dir end

function BasicBc:getEdge() return self._edge end

function BasicBc:label() return "Flux"..self._dirlabel..self._edge end

return BasicBc
