local ffi = require "ffi"

ffi.cdef [[
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};
]]
