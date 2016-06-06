--------------------------------------------------------------------------------
-- Non-uniform Cartesian mesh
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- define C interfaces
ffi.cdef [[
      typedef struct { double _lower, _upper; int32_t _numCells; } Slide_t;
]]

