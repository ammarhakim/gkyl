-- Gkyl ------------------------------------------------------------------------
--
-- Lua interface to gkylzero's gkyl_array structure
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local Alloc = require "Lib.Alloc"
local cuda

_M = {}

ffi.cdef [[
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};

/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
  enum gkyl_elem_type type; // type of data stored in array
  size_t elemsz, ncomp; // size of elements, number of 'components'
  size_t size; // number of indices

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  uint32_t flags;  
  struct gkyl_ref_count ref_count;

  int nthreads, nblocks; // threads per block, number of blocks
  struct gkyl_array *on_dev; // pointer to itself or device data
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with data on NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the
 * host.  However, the on_dev member for this cal is set to a device
 * clone of the host struct, and is what should be used to pass to
 * CUDA kernels which require the entire array struct on device.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Returns true if array lives on NV-GPU.
 *
 * @param arr Array to check
 * @return true of array lives on NV-GPU, false otherwise
 */
bool gkyl_array_is_cu_dev(const struct gkyl_array *arr);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);
]]

local longSz = sizeof("long")

-- Various types for arrays of basic C-types
_M.float = 'GKYL_FLOAT'
_M.double = 'GKYL_DOUBLE'
_M.int = 'GKYL_INT'
_M.long = typeof('long')
_M.char = typeof('char')

local function getArrayTypeCode(atype)
   if atype == typeof("int") then
      return 1
   elseif atype == typeof("float") then
      return 2
   elseif atype == typeof("double") then
      return 3
   end
   return 42 -- user-defined type
end

local function getType(enum)
   if enum == 0 then
      return "int"
   elseif enum == 1 then
      return "float"
   elseif enum == 2 then
      return "double"
   end
end

-- Array ctype
local ArrayCt = typeof("struct gkyl_array")

local array_fn = {
   copy = function (self, src)
      return ffiC.gkyl_array_copy(self, src)
   end,
   clone = function (self)
      return ffiC.gkyl_array_clone(self)
   end,
   aquire = function (self)
      return ffiC.gkyl_array_acquire(self)
   end,
   release = function (self)
      return ffiC.gkyl_array_release(self)
   end,
   fetch = function (self, loc)
      if loc == nil then loc = 0 end
      return ffi.cast(getType(self.type).."*", self.data) + loc*self.ncomp
   end,
   cfetch = function (self, loc)
      if loc == nil then loc = 0 end
      return ffi.cast("const "..getType(self.type).."*", self.data) + loc*self.ncomp
   end,
   get_size = function (self)
      return tonumber(self.size)
   end,
   get_ncomp = function (self)
      return tonumber(self.ncomp)
   end,
   get_elemsz = function (self)
      return tonumber(self.elemsz)
   end,
}

local array_mt = {
   __new = function (self, atype, ncomp, size)
      return ffiC.gkyl_array_new(atype, ncomp, size)
   end,
   __gc = function(self)
      self:release()
   end,
   __index = array_fn,
}
local ArrayCtor = metatype(ArrayCt, array_mt)

-- Construct array of given shape and type
_M.Array = function (atype, ncomp, size)
   return ArrayCtor(atype, ncomp, size)
end

return _M
