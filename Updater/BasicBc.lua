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
local CartDecomp     = require "Lib.CartDecomp"
local Grid           = require "Grid"
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local Proto          = require "Lib.Proto"
local Range          = require "Lib.Range"
local UpdaterBase    = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Constants to represent lower/upper edges
enum gkyl_edge_loc { GKYL_LOWER_EDGE = 0, GKYL_UPPER_EDGE = 1 };

// BC types in this updater.
enum gkyl_bc_basic_type { GKYL_BC_ABSORB = 0, GKYL_BC_REFLECT = 1 };

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

--   if self._grid._zero then
      self._bcType  = assert(tbl.bcType, "Updater.BasicBc: Must specify BC type in 'bcType'.")
      self._basis   = assert(tbl.basis, "Updater.BasicBc: Must specify the basis in 'basis'.")
      local onField = assert(tbl.onField, "Updater.BasicBc: Must specify the field we'll apply BCs to in 'onField'.")

      local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
      local localExtRange = onField:localExtRange()
      local numGhostVec   = self._edge == 'lower' and onField._lowerGhostVec or onField._upperGhostVec

      local bctype -- Match gkyl_bc_basic_type in gkylzero/zero/gkyl_bc_basic.h
      if self._bcType == "absorb" then bctype = 0
      elseif self._bcType == "reflect" then bctype = 1 end

      self._zero = ffi.gc(ffiC.gkyl_bc_basic_new(self._dir-1, edge, localExtRange, numGhostVec:data(), bctype,
                                                 self._basis._zero, onField:numComponents(), cDim, GKYL_USE_GPU or 0),
                          ffiC.gkyl_bc_basic_release)
--   else
--      -- g2 code, to be deleted.
--      self._bcList = assert(
--         tbl.boundaryConditions, "Updater.BasicBc: Must specify boundary conditions to apply with 'boundaryConditions'")
--
--      self._skinLoop = tbl.skinLoop and tbl.skinLoop or "pointwise"
--      if self._skinLoop == "flip" then
--         self._vdir = assert(tbl.vdir, "Updater.BasicBc: Must specify velocity direction to flip with 'vdir'")
--      end
--
--      local advArgs = tbl.advanceArgs  -- Sample arguments for advance method.
--
--      self._idxIn  = Lin.IntVec(self._grid:ndim())
--      self._idxOut = Lin.IntVec(self._grid:ndim())
--      self._xcIn   = Lin.Vec(self._grid:ndim())
--      self._xcOut  = Lin.Vec(self._grid:ndim())
--
--      -- Initialize tools constructed from fields (e.g. ranges).
--      self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil
--   end

   local dirlabel = {"X", "Y", "Z"}
   self._dirlabel = dirlabel[self._dir]
end

function BasicBc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1   -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function BasicBc:initFldTools(inFld, outFld)
   -- Pre-initialize tools (ranges, pointers, etc) depending on fields and used in the advance method.
   local tools = {}

   local qOut  = assert(outFld[1], "BasicBc.advance: Must-specify an output field")

   local grid = self._grid

   local global     = qOut:globalRange()
   local globalExt  = qOut:globalExtRange()
   local localExt   = qOut:localExtRange()
   local ghostRange = localExt:intersect(self:getGhostRange(global, globalExt))   -- Range spanning ghost cells.
   -- Decompose ghost region into threads.
   tools.ghostRangeDecomp = LinearDecomp.LinearDecompRange {
      range = ghostRange, numSplit = grid:numSharedProcs() }

   -- Get the in and out indexers. 
   tools.indexerOut, tools.indexerIn = qOut:genIndexer(), qOut:genIndexer()

   self.flipIdx = self._skinLoop == "flip" 
      and function(idxIn) idxIn[self._vdir] = global:upper(self._vdir) + 1 - idxIn[self._vdir] end
      or function(idxIn) end

   return tools
end

function BasicBc:_advanceBasic(tCurr, inFld, outFld)
   -- Advance method for the basic BCs (copy, absorb, reflect, open, extern).
   local qOut = assert(outFld[1], "BasicBc.advance: Must-specify an output field")

   local grid      = self._grid
   local dir, edge = self._dir, self._edge
   local global    = qOut:globalRange()

   -- Get the in and out pointers.
   local ptrOut, ptrIn = qOut:get(1), qOut:get(1)

   local tId = grid:subGridSharedId() -- Local thread ID.
   for idxOut in self.fldTools.ghostRangeDecomp:rowMajorIter(tId) do 
      qOut:fill(self.fldTools.indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(self._idxIn)
      self.flipIdx(self._idxIn)   -- Reverse the in velocity index if needed
      self._idxIn[dir] = edge=="lower" and global:lower(dir) or global:upper(dir)

      qOut:fill(self.fldTools.indexerIn(self._idxIn), ptrIn)

      self._grid:setIndex(self._idxIn);  self._grid:cellCenter(self._xcIn);
      self._grid:setIndex(idxOut);       self._grid:cellCenter(self._xcOut)

      for _, bc in ipairs(self._bcList) do
         -- Apply the 'bc' function. This can represent many boundary
         -- condition types ranging from a simple copy or a reflection
         -- with the sign flit to QM based electron emission model.
         bc(dir, tCurr, self._idxIn, ptrIn, ptrOut, self._xcOut, self._xcIn)
      end
   end
end

function BasicBc:_advance(tCurr, inFld, outFld)

--   if self._zero then
      local bufferIn = assert(inFld[1], "BasicBc.advance: Must-specify a buffer as large as the ghost cells for this BC.")
      local qOut     = assert(outFld[1], "BasicBc.advance: Must-specify an output field")
      ffiC.gkyl_bc_basic_advance(self._zero, bufferIn._zero, qOut._zero)
--   else
--      -- g2 code, to be deleted.
--      self.fldTools = self.fldTools or self:initFldTools(inFld,outFld)
--
--      self:_advanceBasic(tCurr, inFld, outFld)
--   end
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
