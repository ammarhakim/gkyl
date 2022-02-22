-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local ZeroArray = require "DataStruct.ZeroArray"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local arr = ZeroArray.Array( ZeroArray.double, 1, 200)
   assert_equal(200, arr:get_size(), "Testing array size")
   assert_equal(1, arr:get_ncomp(), "Testing array ncomp")
   assert_equal(8, arr:get_elemsz(), "Testing array elemsz")

   local dat = arr:fetch(0)
   for i=0, arr:get_size()-1 do
      dat[i] = (i+0.5)*0.1
   end

   local brr = arr:clone()
   assert_equal(200, brr:get_size(), "Testing clone array size")
   assert_equal(1, brr:get_ncomp(), "Testing clone array ncomp")
   assert_equal(8, brr:get_elemsz(), "Testing clone array elemsz")

   local dat = brr:fetch(0)
   for i=0, brr:get_size()-1 do
      assert_equal(dat[i], (i+0.5)*0.1)
   end
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
