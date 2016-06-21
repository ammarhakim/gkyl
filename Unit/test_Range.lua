-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local range = Range.Range({0, 0}, {1, 5})

   assert_equal(2, range:ndim(), "Checking dimension")
   assert_equal(12, range:volume(), "Checking volume")   

   assert_equal(0, range:lower(1), "Checking lower")
   assert_equal(0, range:lower(2), "Checking lower")

   assert_equal(1, range:upper(1), "Checking upper")
   assert_equal(5, range:upper(2), "Checking upper")

   assert_equal(2, range:shape(1), "Checking shape")
   assert_equal(6, range:shape(2), "Checking shape")
end

function test_2()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})
   assert_equal(10*20*30, range:volume(), "Checking volume")

   assert_equal(3, range:ndim(), "Checking dimension")

   assert_equal(1, range:lower(1), "Checking lower")
   assert_equal(1, range:lower(2), "Checking lower")
   assert_equal(1, range:lower(3), "Checking lower")   

   assert_equal(10, range:upper(1), "Checking upper")
   assert_equal(20, range:upper(2), "Checking upper")
   assert_equal(30, range:upper(3), "Checking upper")

   assert_equal(10, range:shape(1), "Checking shape")
   assert_equal(20, range:shape(2), "Checking shape")
   assert_equal(30, range:shape(3), "Checking shape")
end

function test_3()
   local range = Range.Range({1}, {10})

   local count = range:lower(1)
   for idx in range:colMajorIter() do
      assert_equal(count, idx[1], "Checking col-major indexer")
      count = count + 1
   end
end

function test_4()
   local range = Range.Range({1, 1}, {5, 6})

   -- fill indices for testing
   local indices = {}
   for j = range:lower(2), range:upper(2) do   
      for i = range:lower(1), range:upper(1) do
	 table.insert(indices, {i,j})
      end
   end

   local count = range:lower(1)
   for idx in range:colMajorIter() do
      assert_equal(indices[count][1], idx[1], "Checking col-major index 1")
      assert_equal(indices[count][2], idx[2], "Checking col-major index 2")
      count = count + 1
   end
end

function test_5()
   local range = Range.Range({1, 1}, {5, 6})

   -- fill indices for testing
   local indices = {}
   for i = range:lower(1), range:upper(1) do   
      for j = range:lower(2), range:upper(2) do
	 table.insert(indices, {i,j})
      end
   end

   local count = range:lower(1)
   for idx in range:rowMajorIter() do
      assert_equal(indices[count][1], idx[1], "Checking col-major index 1")
      assert_equal(indices[count][2], idx[2], "Checking col-major index 2")
      count = count + 1
   end
end

function test_6()
   local range = Range.Range( {1, 1}, {0, 10})
   assert_equal(0, range:volume(), "Checking volume")

   for idx in range:rowMajorIter() do
      print("Should not happen")
   end
end

function test_7()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})
   local extRange = range:extend(1, 2)

   assert_equal(0, extRange:lower(1), "Checking lower extended")
   assert_equal(0, extRange:lower(2), "Checking lower extended")
   assert_equal(0, extRange:lower(3), "Checking lower extended")

   assert_equal(12, extRange:upper(1), "Checking upper extended")
   assert_equal(22, extRange:upper(2), "Checking upper extended")
   assert_equal(32, extRange:upper(3), "Checking upper extended")

   assert_equal(13*23*33, extRange:volume(), "Checking volume extended")
end

function test_8()
   local range = Range.Range({1}, {10})
   local indexer = Range.makeRowMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      assert_equal(count, indexer(i), "Checking indexer")
      count = count + 1
   end
end

function test_9()
   local range = Range.Range({1, 1}, {10, 5})
   local indexer = Range.makeRowMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      for j = range:lower(2), range:upper(2) do
	 assert_equal(count, indexer(i, j), "Checking row-major indexer")
	 count = count + 1
      end
   end
end

function test_10()
   local range = Range.Range({1, 1, -1}, {10, 5, 9})
   local indexer = Range.makeRowMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      for j = range:lower(2), range:upper(2) do
	 for k = range:lower(3), range:upper(3) do
	    assert_equal(count, indexer(i, j, k), "Checking row-major indexer")
	    count = count + 1
	 end
      end
   end
end

function test_11()
   local range = Range.Range({1}, {10})
   local indexer = Range.makeColMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      assert_equal(count, indexer(i), "Checking indexer")
      count = count + 1
   end
end

function test_12()
   local range = Range.Range({1, 1}, {10, 5})
   local indexer = Range.makeColMajorIndexer(range)

   local count = 1
   for j = range:lower(2), range:upper(2) do   
      for i = range:lower(1), range:upper(1) do
	 assert_equal(count, indexer(i, j), "Checking col-major indexer")
	 count = count + 1
      end
   end
end

function test_13()
   local range = Range.Range({1, 1, -1}, {10, 5, 9})
   local indexer = Range.makeColMajorIndexer(range)

   local count = 1
   for k = range:lower(3), range:upper(3) do
      for j = range:lower(2), range:upper(2) do      
	 for i = range:lower(1), range:upper(1) do
	    assert_equal(count, indexer(i, j, k), "Checking col-major indexer")
	    count = count + 1
	 end
      end
   end
end

-- Run tests
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
test_7()
test_8()
test_9()
test_10()
test_11()
test_12()
test_13()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
