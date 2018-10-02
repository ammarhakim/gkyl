-- Gkyl ------------------------------------------------------------------------
--
-- Test for linear decomposition object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local LinearDecomp = require "Lib.LinearDecomp"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_0()
   local linDecomp = LinearDecomp.LinearDecomp { domSize = 100, numSplit = 1 }

   assert_equal(100, linDecomp:domSize(), "Checking dom size")
   assert_equal(1, linDecomp:numSplit(), "Checking num splits")

   local count = 0
   for d = 1, linDecomp:numSplit() do
      count = count+linDecomp:shape(d)
   end
   assert_equal(linDecomp:domSize(), count, "Checking if total shape is correct")

   for d = 1, linDecomp:numSplit() do
      assert_equal(100, linDecomp:shape(d), "Checking local shape")
      assert_equal(100*(d-1)+1, linDecomp:lower(d), "Checking local start")
      assert_equal(100*d, linDecomp:upper(d), "Checking local end")
   end   
end

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

function test_3()
   local r = Range.Range({1}, {100})
   local linDecomp = LinearDecomp.LinearDecompRange { range = r, numSplit = 10 }

   for d = 1, linDecomp:numSplit() do
      assert_equal(10, linDecomp:shape(d), "Checking size")
      assert_equal((d-1)*10+1, linDecomp:lower(d), "Checking lower")
      assert_equal((d-1)*10+10, linDecomp:upper(d), "Checking upper")
      
      local idx = linDecomp:rowStartIndex(d)
      assert_equal((d-1)*10+1, idx[1], "Checking index")
   end

   for d = 1, linDecomp:numSplit() do
      local idx = linDecomp:colStartIndex(d)
      assert_equal((d-1)*10+1, idx[1], "Checking index")
   end   
end

function test_4()
   local r = Range.Range({1, 1}, {10, 10})
   local linDecomp = LinearDecomp.LinearDecompRange { range = r, numSplit = 10 }

   local count = 0
   for d = 1, linDecomp:numSplit() do
      local startIdx, sz = linDecomp:colStartIndex(d), linDecomp:shape(d)
      for idx in r:colMajorIter(startIdx, sz) do
	 count = count+1
      end
   end
   assert_equal(r:volume(), count, "Checking volume")

   local count = 0
   for d = 1, linDecomp:numSplit() do
      local startIdx, sz = linDecomp:rowStartIndex(d), linDecomp:shape(d)
      for idx in r:rowMajorIter(startIdx, sz) do
	 count = count+1
      end
   end
   assert_equal(r:volume(), count, "Checking volume")
end

function test_5()
   local linDecomp = LinearDecomp.LinearDecomp { domSize = 16, numSplit = 32 }

   for d = 1, linDecomp:domSize() do
      assert_equal(1, linDecomp:shape(d), "Checking splitting 16 into 32")
   end
   for d = linDecomp:domSize()+1, linDecomp:numSplit() do
      assert_equal(0, linDecomp:shape(d), "Checking splitting 16 into 32")
   end   
end

function test_6()
   local r = Range.Range({1, 1}, {10, 10})
   local linDecomp = LinearDecomp.LinearDecompRange { range = r, numSplit = 10 }

   local count = 0
   for t = 1, linDecomp:numSplit() do
      for itr in linDecomp:colMajorIter(t) do
	 count = count+1
      end
   end
   assert_equal(r:volume(), count, "Checking volume")

   count = 0
   for t = 1, linDecomp:numSplit() do
      for itr in linDecomp:rowMajorIter(t) do
	 count = count+1
      end
   end
   assert_equal(r:volume(), count, "Checking volume")   
end

-- Run tests
test_0()
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
