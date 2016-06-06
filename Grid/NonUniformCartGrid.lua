--------------------------------------------------------------------------------
-- Non-uniform Cartesian mesh
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- define C interfaces
ffi.cdef [[
      typedef struct { double _lower, _upper; int32_t _numCells; } Slide_t;
]]

