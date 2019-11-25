-- Gkyl ------------------------------------------------------------------------
--
-- Test for memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Lin = require "Lib.Linalg"
local Unit = require "Unit"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
]]

function test_1()
   local x, y = Lin.Vec(10), Lin.Vec(10)
   for i = 1, #y do
      y[i] = 1.5+i
      x[i] = 2.5+10*i
   end

   ffi.C.cblas_daxpy(#y, 1.0, x:data(), 1, y:data(), 1)

   for i = 1, #y do
      assert_equal(1.5+i+2.5+10*i, y[i], "Checking daxypy")
   end
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
