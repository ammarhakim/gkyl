-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"

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

   local extRange1 = range:extendDir(1, 1, 2)
   assert_equal(0, extRange1:lower(1), "Checking lower extended")
   assert_equal(1, extRange1:lower(2), "Checking lower extended")
   assert_equal(1, extRange1:lower(3), "Checking lower extended")

   assert_equal(12, extRange1:upper(1), "Checking upper extended")
   assert_equal(20, extRange1:upper(2), "Checking upper extended")
   assert_equal(30, extRange1:upper(3), "Checking upper extended")

   extRange1 = range:extendDir(2, 1, 2)
   assert_equal(1, extRange1:lower(1), "Checking lower extended")
   assert_equal(0, extRange1:lower(2), "Checking lower extended")
   assert_equal(1, extRange1:lower(3), "Checking lower extended")

   assert_equal(10, extRange1:upper(1), "Checking upper extended")
   assert_equal(22, extRange1:upper(2), "Checking upper extended")
   assert_equal(30, extRange1:upper(3), "Checking upper extended")

   extRange1 = range:extendDir(3, 1, 2)
   assert_equal(1, extRange1:lower(1), "Checking lower extended")
   assert_equal(1, extRange1:lower(2), "Checking lower extended")
   assert_equal(0, extRange1:lower(3), "Checking lower extended")

   assert_equal(10, extRange1:upper(1), "Checking upper extended")
   assert_equal(20, extRange1:upper(2), "Checking upper extended")
   assert_equal(32, extRange1:upper(3), "Checking upper extended")
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

function test_14()
   local r1 = Range.Range({1, 1}, {10, 10})
   local r2 = Range.Range({1, 1}, {12, 10})
   local r3 = Range.Range({1}, {10})
   local r4 = Range.Range({1, 1}, {10, 10})
   
   assert_equal(false, r1 == r2, "Checking for range equality")
   assert_equal(false, r1 == r3, "Checking for range equality")
   assert_equal(true, r1 == r4, "Checking for range equality")

   assert_equal(true, r1 ~= r2, "Checking for range inequality")
   assert_equal(true, r1 ~= r3, "Checking for range inequality")
   assert_equal(false, r1 ~= r4, "Checking for range inequality")

end

function test_15()
   local bigr = Range.Range({1, 1}, {10, 10})
   local range = bigr:shorten(1)

   assert_equal(1, range:lower(1), "Checking shorter range")
   assert_equal(1, range:upper(1), "Checking shorter range")
   assert_equal(1, range:shape(1), "Checking shorter range")   

   assert_equal(1, range:lower(2), "Checking shorter range")
   assert_equal(10, range:upper(2), "Checking shorter range")
   assert_equal(10, range:shape(2), "Checking shorter range")

   assert_equal(10*1, range:volume(), "Checking shorter volume")
end

function test_16()
   local r = Range.Range({1, 2, 3}, {10, 11, 12})
   local lv, uv = r:lowerAsVec(), r:upperAsVec()
   for i = 1, 3 do
      assert_equal(i, lv[i], "Checking lower vector")
      assert_equal(9+i, uv[i], "Checking upper vector")
   end
end

function test_17()
   local r1 = Range.Range({1, 2}, {10, 11})

   local r2 = Range.Range({5, 6}, {12, 13})
   local r3 = r1:intersect(r2)
   assert_equal(5, r3:lower(1), "Checking intersected range")
   assert_equal(6, r3:lower(2), "Checking intersected range")
   assert_equal(10, r3:upper(1), "Checking intersected range")
   assert_equal(11, r3:upper(2), "Checking intersected range")

   local r3 = r2:intersect(r1)
   assert_equal(5, r3:lower(1), "Checking intersected range")
   assert_equal(6, r3:lower(2), "Checking intersected range")
   assert_equal(10, r3:upper(1), "Checking intersected range")
   assert_equal(11, r3:upper(2), "Checking intersected range")

   local r2 = Range.Range({5, 5}, {10, 10})
   local r3 = r1:intersect(r2)
   assert_equal(5, r3:lower(1), "Checking intersected range")
   assert_equal(5, r3:lower(2), "Checking intersected range")
   assert_equal(10, r3:upper(1), "Checking intersected range")
   assert_equal(10, r3:upper(2), "Checking intersected range")

   local r2 = Range.Range({4, 4}, {15, 8})
   local r3 = r1:intersect(r2)
   assert_equal(4, r3:lower(1), "Checking intersected range")
   assert_equal(4, r3:lower(2), "Checking intersected range")
   assert_equal(10, r3:upper(1), "Checking intersected range")
   assert_equal(8, r3:upper(2), "Checking intersected range")

   local r2 = Range.Range({20, 20}, {25, 25})
   local r3 = r1:intersect(r2)
   assert_equal(0, r3:volume(), "Checking volume of empty intersection")

end

function test_18()
   local r1 = Range.Range({1, 2}, {10, 11})
   local offsets = {10, 12}
   local rshift = r1:shift(offsets)

   assert_equal(r1:volume(), rshift:volume(), "Checking shifted volumes")
   for d = 1, r1:ndim() do
      assert_equal(rshift:lower(d)-r1:lower(d), offsets[d], "Checking offset")
      assert_equal(rshift:upper(d)-r1:upper(d), offsets[d], "Checking offset")
   end

   offsets = {-10, -12}
   rshift = r1:shift(offsets)

   assert_equal(r1:volume(), rshift:volume(), "Checking shifted volumes")
   for d = 1, r1:ndim() do
      assert_equal(rshift:lower(d)-r1:lower(d), offsets[d], "Checking offset")
      assert_equal(rshift:upper(d)-r1:upper(d), offsets[d], "Checking offset")
   end
end

function test_19()
   local r1 = Range.Range({1, 2}, {10, 11})
   local offset = 10
   local rshift = r1:shiftInDir(1, offset)

   assert_equal(r1:volume(), rshift:volume(), "Checking shifted volumes")
   assert_equal(rshift:lower(1)-r1:lower(1), offset, "Checking offset")
   assert_equal(rshift:upper(1)-r1:upper(1), offset, "Checking offset")

   offset = -12
   rshift = r1:shiftInDir(2, offset)

   assert_equal(r1:volume(), rshift:volume(), "Checking shifted volumes")
   assert_equal(rshift:lower(2)-r1:lower(2), offset, "Checking offset")
   assert_equal(rshift:upper(2)-r1:upper(2), offset, "Checking offset")
end

function test_20()
   local r = Range.Range({1, 1}, {10, 10})
   local count = 0
   for idx in r:colMajorIter() do
      count = count+1
   end
   assert_equal(r:volume(), count, "Checking if iterator bumped over full range")

   count = 0
   for idx in r:rowMajorIter() do
      count = count+1
   end
   assert_equal(r:volume(), count, "Checking if iterator bumped over full range")   
end

function test_21()
   local r = Range.Range({1, 1}, {10, 10})
   local sidx = Lin.IntVec(r:ndim())

   -- defaults
   for d = 1, r:ndim() do
      sidx[d] = r:lower(d)
   end

   local count = 0
   for idx in r:colMajorIter(sidx, r:volume()) do
      count = count+1
   end
   assert_equal(r:volume(), count, "Checking if iterator bumped over full range")

   -- starting at (5,5) over remaining domain
   sidx[1] = 5; sidx[2] = 5
   count = 0
   for idx in r:rowMajorIter(sidx) do
      count = count+1
   end
   assert_equal(56, count, "Checking if iterator bumped over full range")

   count = 0
   for idx in r:colMajorIter(sidx) do
      count = count+1
   end
   assert_equal(56, count, "Checking if iterator bumped over full range")

   -- starting at (5,5) over 10 cells
   count = 0
   for idx in r:rowMajorIter(sidx, 10) do
      count = count+1
   end
   assert_equal(10, count, "Checking if iterator bumped over full range")

   -- starting at (5,5) over 0 cells
   count = 0
   for idx in r:rowMajorIter(sidx, 0) do
      count = count+1
   end
   assert_equal(0, count, "Checking if iterator bumped over full range")
end

function test_22()
   local r = Range.Range({-1}, {10})
   local invIndexer = Range.makeRowMajorInvIndexer(r)
   local idx = Lin.IntVec(r:ndim())

   local ix = r:lower(1)
   for loc = 1, r:volume() do
      invIndexer(loc, idx)
      assert_equal(ix, idx[1], "Checking inverse indexer")
      ix = ix+1
   end
end

function test_23()
   local r = Range.Range({-1, -5}, {1, 10})
   local indexer = Range.makeRowMajorGenIndexer(r)
   local invIndexer = Range.makeRowMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:rowMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "Checking inv indexer")
   end
end

function test_24()
   local r = Range.Range({-1, -5, 0}, {1, 10, 4})
   local indexer = Range.makeRowMajorGenIndexer(r)
   local invIndexer = Range.makeRowMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:rowMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "Checking inv indexer")
      assert_equal(idx[3], invIdx[3], "Checking inv indexer")
   end
end

function test_25()
   local r = Range.Range({-1}, {10})
   local invIndexer = Range.makeColMajorInvIndexer(r)
   local idx = Lin.IntVec(r:ndim())

   local ix = r:lower(1)
   for loc = 1, r:volume() do
      invIndexer(loc, idx)
      assert_equal(ix, idx[1], "Checking inverse indexer")
      ix = ix+1
   end
end

function test_26()
   local r = Range.Range({-1, -5}, {1, 10})
   local indexer = Range.makeColMajorGenIndexer(r)
   local invIndexer = Range.makeColMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:colMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "Checking inv indexer")
   end
end

function test_27()
   local r = Range.Range({-1, -5, 0}, {1, 10, 4})
   local indexer = Range.makeColMajorGenIndexer(r)
   local invIndexer = Range.makeColMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:colMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "Checking inv indexer")
      assert_equal(idx[3], invIdx[3], "Checking inv indexer")
   end
end

function test_28()
   local r = Range.Range({1, 1}, {10, 10})

   local lowerSkin = r:lowerSkin(1, 1)
   assert_equal(10, lowerSkin:volume(), "Checking volume")
   assert_equal(1, lowerSkin:lower(1), "Checking lower in 1")
   assert_equal(1, lowerSkin:upper(1), "Checking upper in 1")
   assert_equal(1, lowerSkin:lower(2), "Checking lower in 2")
   assert_equal(10, lowerSkin:upper(2), "Checking upper in 2")

   local lowerSkin = r:lowerSkin(2, 1)
   assert_equal(10, lowerSkin:volume(), "Checking volume")
   assert_equal(1, lowerSkin:lower(1), "Checking lower in 1")
   assert_equal(10, lowerSkin:upper(1), "Checking upper in 1")
   assert_equal(1, lowerSkin:lower(2), "Checking lower in 2")
   assert_equal(1, lowerSkin:upper(2), "Checking upper in 2")

   local upperSkin = r:upperSkin(1, 1)
   assert_equal(10, upperSkin:volume(), "Checking volume")
   assert_equal(10, upperSkin:lower(1), "Checking upper in 1")
   assert_equal(10, upperSkin:upper(1), "Checking upper in 1")
   assert_equal(1, upperSkin:lower(2), "Checking upper in 2")
   assert_equal(10, upperSkin:upper(2), "Checking upper in 2")

   local upperSkin = r:upperSkin(2, 1)
   assert_equal(10, upperSkin:volume(), "Checking volume")
   assert_equal(1, upperSkin:lower(1), "Checking upper in 1")
   assert_equal(10, upperSkin:upper(1), "Checking upper in 1")
   assert_equal(10, upperSkin:lower(2), "Checking upper in 2")
   assert_equal(10, upperSkin:upper(2), "Checking upper in 2")

   ---

   local lowerGhost = r:lowerGhost(1, 1)
   assert_equal(10, lowerGhost:volume(), "Checking volume")
   assert_equal(0, lowerGhost:lower(1), "Checking lower in 1")
   assert_equal(0, lowerGhost:upper(1), "Checking upper in 1")
   assert_equal(1, lowerGhost:lower(2), "Checking lower in 2")
   assert_equal(10, lowerGhost:upper(2), "Checking upper in 2")

   local lowerGhost = r:lowerGhost(2, 1)
   assert_equal(10, lowerGhost:volume(), "Checking volume")
   assert_equal(1, lowerGhost:lower(1), "Checking lower in 1")
   assert_equal(10, lowerGhost:upper(1), "Checking upper in 1")
   assert_equal(0, lowerGhost:lower(2), "Checking lower in 2")
   assert_equal(0, lowerGhost:upper(2), "Checking upper in 2")

   local upperGhost = r:upperGhost(1, 1)
   assert_equal(10, upperGhost:volume(), "Checking volume")
   assert_equal(11, upperGhost:lower(1), "Checking upper in 1")
   assert_equal(11, upperGhost:upper(1), "Checking upper in 1")
   assert_equal(1, upperGhost:lower(2), "Checking upper in 2")
   assert_equal(10, upperGhost:upper(2), "Checking upper in 2")

   local upperGhost = r:upperGhost(2, 1)
   assert_equal(10, upperGhost:volume(), "Checking volume")
   assert_equal(1, upperGhost:lower(1), "Checking upper in 1")
   assert_equal(10, upperGhost:upper(1), "Checking upper in 1")
   assert_equal(11, upperGhost:lower(2), "Checking upper in 2")
   assert_equal(11, upperGhost:upper(2), "Checking upper in 2")

end

function test_29()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})

   local extRange1 = range:extendDirs({1}, 1, 2)
   assert_equal(0, extRange1:lower(1), "Checking lower extended")
   assert_equal(1, extRange1:lower(2), "Checking lower extended")
   assert_equal(1, extRange1:lower(3), "Checking lower extended")

   assert_equal(12, extRange1:upper(1), "Checking upper extended")
   assert_equal(20, extRange1:upper(2), "Checking upper extended")
   assert_equal(30, extRange1:upper(3), "Checking upper extended")

   extRange1 = range:extendDirs({2}, 1, 2)
   assert_equal(1, extRange1:lower(1), "Checking lower extended")
   assert_equal(0, extRange1:lower(2), "Checking lower extended")
   assert_equal(1, extRange1:lower(3), "Checking lower extended")

   assert_equal(10, extRange1:upper(1), "Checking upper extended")
   assert_equal(22, extRange1:upper(2), "Checking upper extended")
   assert_equal(30, extRange1:upper(3), "Checking upper extended")

   extRange1 = range:extendDirs({3}, 1, 2)
   assert_equal(1, extRange1:lower(1), "Checking lower extended")
   assert_equal(1, extRange1:lower(2), "Checking lower extended")
   assert_equal(0, extRange1:lower(3), "Checking lower extended")

   assert_equal(10, extRange1:upper(1), "Checking upper extended")
   assert_equal(20, extRange1:upper(2), "Checking upper extended")
   assert_equal(32, extRange1:upper(3), "Checking upper extended")

   extRange1 = range:extendDirs({1,3}, 1, 2)
   assert_equal(0, extRange1:lower(1), "Checking lower extended")
   assert_equal(1, extRange1:lower(2), "Checking lower extended")
   assert_equal(0, extRange1:lower(3), "Checking lower extended")

   assert_equal(12, extRange1:upper(1), "Checking upper extended")
   assert_equal(20, extRange1:upper(2), "Checking upper extended")
   assert_equal(32, extRange1:upper(3), "Checking upper extended")   
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
test_14()
test_15()
test_16()
test_17()
test_18()
test_19()
test_20()
test_21()
test_22()
test_23()
test_24()
test_25()
test_26()
test_27()
test_28()
test_29()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
