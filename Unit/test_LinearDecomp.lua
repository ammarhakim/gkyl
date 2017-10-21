-- Gkyl ------------------------------------------------------------------------
--
-- Test for linear decomposition object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local LinearDecomp = require "Lib.LinearDecomp"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local linDecomp = LinearDecomp.LinearDecomp { domSize = 100, numSplit = 10 }

   assert_equal(100, linDecomp:domSize(), "Checking dom size")
   assert_equal(10, linDecomp:numSplit(), "Checking num splits")

   local count = 0
   for d = 1, linDecomp:numSplit() do
      count = count+linDecomp:shape(d)
   end
   assert_equal(linDecomp:domSize(), count, "Checking if total shape is correct")

   for d = 1, linDecomp:numSplit() do
      assert_equal(10, linDecomp:shape(d), "Checking local shape")
      assert_equal(10*(d-1)+1, linDecomp:lower(d), "Checking local start")
      assert_equal(10*d, linDecomp:upper(d), "Checking local end")
   end   
end

function test_2()
   local linDecomp = LinearDecomp.LinearDecomp { domSize = 100, numSplit = 7 }

   assert_equal(100, linDecomp:domSize(), "Checking dom size")
   assert_equal(7, linDecomp:numSplit(), "Checking num splits")

   local count = 0
   for d = 1, linDecomp:numSplit() do
      count = count+linDecomp:shape(d)
   end
   assert_equal(linDecomp:domSize(), count, "Checking if total shape is correct")
end

-- Run tests
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
