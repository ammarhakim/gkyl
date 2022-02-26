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

ffi.cdef [[
  void unit_showRange(struct gkyl_range *range);
]]

function test_1()
   local range = Range.Range({0, 0}, {1, 5})

   assert_equal(2, range:ndim(), "test_1: Checking dimension")
   assert_equal(12, range:volume(), "test_1: Checking volume")   

   assert_equal(0, range:lower(1), "test_1: Checking lower")
   assert_equal(0, range:lower(2), "test_1: Checking lower")

   assert_equal(1, range:upper(1), "test_1: Checking upper")
   assert_equal(5, range:upper(2), "test_1: Checking upper")

   assert_equal(2, range:shape(1), "test_1: Checking shape")
   assert_equal(6, range:shape(2), "test_1: Checking shape")

end

function test_2()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})
   assert_equal(10*20*30, range:volume(), "test_1: Checking volume")

   assert_equal(3, range:ndim(), "test_1: Checking dimension")

   assert_equal(1, range:lower(1), "test_1: Checking lower")
   assert_equal(1, range:lower(2), "test_1: Checking lower")
   assert_equal(1, range:lower(3), "test_1: Checking lower")   

   assert_equal(10, range:upper(1), "test_1: Checking upper")
   assert_equal(20, range:upper(2), "test_1: Checking upper")
   assert_equal(30, range:upper(3), "test_1: Checking upper")

   assert_equal(10, range:shape(1), "test_1: Checking shape")
   assert_equal(20, range:shape(2), "test_1: Checking shape")
   assert_equal(30, range:shape(3), "test_1: Checking shape")
end

function test_3()
   local range = Range.Range({1}, {10})

   local count = range:lower(1)
   for idx in range:colMajorIter() do
      assert_equal(count, idx[1], "test_3: Checking col-major indexer")
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
      assert_equal(indices[count][1], idx[1], "test_4: Checking col-major index 1")
      assert_equal(indices[count][2], idx[2], "test_4: Checking col-major index 2")
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
      assert_equal(indices[count][1], idx[1], "test_5: Checking col-major index 1")
      assert_equal(indices[count][2], idx[2], "test_5: Checking col-major index 2")
      count = count + 1
   end
end

function test_6()
   local range = Range.Range( {1, 1}, {0, 10})
   assert_equal(0, range:volume(), "test_6: Checking volume")

   for idx in range:rowMajorIter() do
      print("Should not happen")
   end
end

function test_7()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})

   local rangeCopy = Range.Range({0},{1})
   rangeCopy:copy(range)

   assert_equal(1, rangeCopy:lower(1), "test_7: Checking lower copy")
   assert_equal(1, rangeCopy:lower(2), "test_7: Checking lower copy")
   assert_equal(1, rangeCopy:lower(3), "test_7: Checking lower copy")

   assert_equal(10, rangeCopy:upper(1), "test_7: Checking upper copy")
   assert_equal(20, rangeCopy:upper(2), "test_7: Checking upper copy")
   assert_equal(30, rangeCopy:upper(3), "test_7: Checking upper copy")

   local extRange = range:extend(1, 2)

   assert_equal(0, extRange:lower(1), "test_7: Checking lower extended")
   assert_equal(0, extRange:lower(2), "test_7: Checking lower extended")
   assert_equal(0, extRange:lower(3), "test_7: Checking lower extended")

   assert_equal(12, extRange:upper(1), "test_7: Checking upper extended")
   assert_equal(22, extRange:upper(2), "test_7: Checking upper extended")
   assert_equal(32, extRange:upper(3), "test_7: Checking upper extended")

   assert_equal(13*23*33, extRange:volume(), "test_7: Checking volume extended")

   local extRange1 = range:extendDir(1, 1, 2)
   assert_equal(0, extRange1:lower(1), "test_7: Checking lower extended")
   assert_equal(1, extRange1:lower(2), "test_7: Checking lower extended")
   assert_equal(1, extRange1:lower(3), "test_7: Checking lower extended")

   assert_equal(12, extRange1:upper(1), "test_7: Checking upper extended")
   assert_equal(20, extRange1:upper(2), "test_7: Checking upper extended")
   assert_equal(30, extRange1:upper(3), "test_7: Checking upper extended")

   extRange1 = range:extendDir(2, 1, 2)
   assert_equal(1, extRange1:lower(1), "test_7: Checking lower extended")
   assert_equal(0, extRange1:lower(2), "test_7: Checking lower extended")
   assert_equal(1, extRange1:lower(3), "test_7: Checking lower extended")

   assert_equal(10, extRange1:upper(1), "test_7: Checking upper extended")
   assert_equal(22, extRange1:upper(2), "test_7: Checking upper extended")
   assert_equal(30, extRange1:upper(3), "test_7: Checking upper extended")

   extRange1 = range:extendDir(3, 1, 2)
   assert_equal(1, extRange1:lower(1), "test_7: Checking lower extended")
   assert_equal(1, extRange1:lower(2), "test_7: Checking lower extended")
   assert_equal(0, extRange1:lower(3), "test_7: Checking lower extended")

   assert_equal(10, extRange1:upper(1), "test_7: Checking upper extended")
   assert_equal(20, extRange1:upper(2), "test_7: Checking upper extended")
   assert_equal(32, extRange1:upper(3), "test_7: Checking upper extended")
end

function test_8()
   local range = Range.Range({1}, {10})
   local indexer = Range.makeRowMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      assert_equal(count, indexer(i), "test_8: Checking indexer")
      count = count + 1
   end
end

function test_9()
   local range = Range.Range({1, 1}, {10, 5})
   local indexer = Range.makeRowMajorIndexer(range)

   local count = 1
   for i = range:lower(1), range:upper(1) do
      for j = range:lower(2), range:upper(2) do
	 assert_equal(count, indexer(i, j), "test_9: Checking row-major indexer")
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
	    assert_equal(count, indexer(i, j, k), "test_10: Checking row-major indexer")
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
      assert_equal(count, indexer(i), "test_11: Checking indexer")
      count = count + 1
   end
end

function test_12()
   local range = Range.Range({1, 1}, {10, 5})
   local indexer = Range.makeColMajorIndexer(range)

   local count = 1
   for j = range:lower(2), range:upper(2) do   
      for i = range:lower(1), range:upper(1) do
	 assert_equal(count, indexer(i, j), "test_12: Checking col-major indexer")
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
	    assert_equal(count, indexer(i, j, k), "test_13: Checking col-major indexer")
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
   
   assert_equal(false, r1 == r2, "test_14: Checking for range equality")
   assert_equal(false, r1 == r3, "test_14: Checking for range equality")
   assert_equal(true, r1 == r4, "test_14: Checking for range equality")

   assert_equal(true, r1 ~= r2, "test_14: Checking for range inequality")
   assert_equal(true, r1 ~= r3, "test_14: Checking for range inequality")
   assert_equal(false, r1 ~= r4, "test_14: Checking for range inequality")

end

function test_15()
   local bigr = Range.Range({1, 1}, {10, 10})
   local range = bigr:shorten(1)

   assert_equal(1, range:lower(1), "test_15: Checking shorter range")
   assert_equal(1, range:upper(1), "test_15: Checking shorter range")
   assert_equal(1, range:shape(1), "test_15: Checking shorter range")   

   assert_equal(1, range:lower(2), "test_15: Checking shorter range")
   assert_equal(10, range:upper(2), "test_15: Checking shorter range")
   assert_equal(10, range:shape(2), "test_15: Checking shorter range")

   assert_equal(10*1, range:volume(), "test_15: Checking shorter volume")

   local range = bigr:shortenFromBelow(1)

   assert_equal(10, range:lower(1), "test_15: Checking shorter range")
   assert_equal(10, range:upper(1), "test_15: Checking shorter range")
   assert_equal(1, range:shape(1), "test_15: Checking shorter range")   

   assert_equal(1, range:lower(2), "test_15: Checking shorter range")
   assert_equal(10, range:upper(2), "test_15: Checking shorter range")
   assert_equal(10, range:shape(2), "test_15: Checking shorter range")

   assert_equal(10*1, range:volume(), "test_15: Checking shorter volume")
end

function test_15_b()
   local bigr = Range.Range({1, 1}, {10, 10})
   local range = bigr:shorten(1, 2)

   assert_equal(1, range:lower(1), "test_15_b: Checking shorter range")
   assert_equal(2, range:upper(1), "test_15_b: Checking shorter range")
   assert_equal(2, range:shape(1), "test_15_b: Checking shorter range")

   assert_equal(1, range:lower(2), "test_15_b: Checking shorter range")
   assert_equal(10, range:upper(2), "test_15_b: Checking shorter range")
   assert_equal(10, range:shape(2), "test_15_b: Checking shorter range")

   assert_equal(10*2, range:volume(), "test_15_b: Checking shorter volume")

   local range = bigr:shortenFromBelow(1, 2)

   assert_equal(9, range:lower(1), "test_15_b: Checking shorter range")
   assert_equal(10, range:upper(1), "test_15_b: Checking shorter range")
   assert_equal(2, range:shape(1), "test_15_b: Checking shorter range")

   assert_equal(1, range:lower(2), "test_15_b: Checking shorter range")
   assert_equal(10, range:upper(2), "test_15_b: Checking shorter range")
   assert_equal(10, range:shape(2), "test_15_b: Checking shorter range")

   assert_equal(10*2, range:volume(), "test_15_b: Checking shorter volume")
end

function test_16()
   local r = Range.Range({1, 2, 3}, {10, 11, 12})
   local lv, uv = r:lowerAsVec(), r:upperAsVec()
   for i = 1, 3 do
      assert_equal(i, lv[i], "test_16: Checking lower vector")
      assert_equal(9+i, uv[i], "test_16: Checking upper vector")
   end
end

function test_17()
   local r1 = Range.Range({1, 2}, {10, 11})

   local r2 = Range.Range({5, 6}, {12, 13})
   local r3 = r1:intersect(r2)
   assert_equal(5, r3:lower(1), "test_17: Checking intersected range")
   assert_equal(6, r3:lower(2), "test_17: Checking intersected range")
   assert_equal(10, r3:upper(1), "test_17: Checking intersected range")
   assert_equal(11, r3:upper(2), "test_17: Checking intersected range")

   local r3 = r2:intersect(r1)
   assert_equal(5, r3:lower(1), "test_17: Checking intersected range")
   assert_equal(6, r3:lower(2), "test_17: Checking intersected range")
   assert_equal(10, r3:upper(1), "test_17: Checking intersected range")
   assert_equal(11, r3:upper(2), "test_17: Checking intersected range")

   local r2 = Range.Range({5, 5}, {10, 10})
   local r3 = r1:intersect(r2)
   assert_equal(5, r3:lower(1), "test_17: Checking intersected range")
   assert_equal(5, r3:lower(2), "test_17: Checking intersected range")
   assert_equal(10, r3:upper(1), "test_17: Checking intersected range")
   assert_equal(10, r3:upper(2), "test_17: Checking intersected range")

   local r2 = Range.Range({4, 4}, {15, 8})
   local r3 = r1:intersect(r2)
   assert_equal(4, r3:lower(1), "test_17: Checking intersected range")
   assert_equal(4, r3:lower(2), "test_17: Checking intersected range")
   assert_equal(10, r3:upper(1), "test_17: Checking intersected range")
   assert_equal(8, r3:upper(2), "test_17: Checking intersected range")

   local r2 = Range.Range({20, 20}, {25, 25})
   local r3 = r1:intersect(r2)
   assert_equal(0, r3:volume(), "test_17: Checking volume of empty intersection")

end

function test_18()
   local r1 = Range.Range({1, 2}, {10, 11})
   local offsets = {10, 12}
   local rshift = r1:shift(offsets)

   assert_equal(r1:volume(), rshift:volume(), "test_18: Checking shifted volumes")
   for d = 1, r1:ndim() do
      assert_equal(rshift:lower(d)-r1:lower(d), offsets[d], "test_18: Checking offset")
      assert_equal(rshift:upper(d)-r1:upper(d), offsets[d], "test_18: Checking offset")
   end

   offsets = {-10, -12}
   rshift = r1:shift(offsets)

   assert_equal(r1:volume(), rshift:volume(), "test_18: Checking shifted volumes")
   for d = 1, r1:ndim() do
      assert_equal(rshift:lower(d)-r1:lower(d), offsets[d], "test_18: Checking offset")
      assert_equal(rshift:upper(d)-r1:upper(d), offsets[d], "test_18: Checking offset")
   end
end

function test_19()
   local r1 = Range.Range({1, 2}, {10, 11})
   local offset = 10
   local rshift = r1:shiftInDir(1, offset)

   assert_equal(r1:volume(), rshift:volume(), "test_19: Checking shifted volumes")
   assert_equal(rshift:lower(1)-r1:lower(1), offset, "test_19: Checking offset")
   assert_equal(rshift:upper(1)-r1:upper(1), offset, "test_19: Checking offset")

   offset = -12
   rshift = r1:shiftInDir(2, offset)

   assert_equal(r1:volume(), rshift:volume(), "test_19: Checking shifted volumes")
   assert_equal(rshift:lower(2)-r1:lower(2), offset, "test_19: Checking offset")
   assert_equal(rshift:upper(2)-r1:upper(2), offset, "test_19: Checking offset")
end

function test_20()
   local r = Range.Range({1, 1}, {10, 10})
   local count = 0
   for idx in r:colMajorIter() do
      count = count+1
   end
   assert_equal(r:volume(), count, "test_20: Checking if iterator bumped over full range")

   count = 0
   for idx in r:rowMajorIter() do
      count = count+1
   end
   assert_equal(r:volume(), count, "test_20: Checking if iterator bumped over full range")   
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
   assert_equal(r:volume(), count, "test_21: Checking if iterator bumped over full range")

   -- starting at (5,5) over remaining domain
   sidx[1] = 5; sidx[2] = 5
   count = 0
   for idx in r:rowMajorIter(sidx) do
      count = count+1
   end
   assert_equal(56, count, "test_21: Checking if iterator bumped over full range")

   count = 0
   for idx in r:colMajorIter(sidx) do
      count = count+1
   end
   assert_equal(56, count, "test_21: Checking if iterator bumped over full range")

   -- starting at (5,5) over 10 cells
   count = 0
   for idx in r:rowMajorIter(sidx, 10) do
      count = count+1
   end
   assert_equal(10, count, "test_21: Checking if iterator bumped over full range")

   -- starting at (5,5) over 0 cells
   count = 0
   for idx in r:rowMajorIter(sidx, 0) do
      count = count+1
   end
   assert_equal(0, count, "test_21: Checking if iterator bumped over full range")
end

function test_22()
   local r = Range.Range({-1}, {10})
   local invIndexer = Range.makeRowMajorInvIndexer(r)
   local idx = Lin.IntVec(r:ndim())

   local ix = r:lower(1)
   for loc = 1, r:volume() do
      invIndexer(loc, idx)
      assert_equal(ix, idx[1], "test_22: Checking inverse indexer")
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
      assert_equal(idx[1], invIdx[1], "test_23: Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "test_23: Checking inv indexer")
   end
end

function test_24()
   local r = Range.Range({-1, -5, 0}, {1, 10, 4})
   local indexer = Range.makeRowMajorGenIndexer(r)
   local invIndexer = Range.makeRowMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:rowMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "test_24: Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "test_24: Checking inv indexer")
      assert_equal(idx[3], invIdx[3], "test_24: Checking inv indexer")
   end
end

function test_25()
   local r = Range.Range({-1}, {10})
   local invIndexer = Range.makeColMajorInvIndexer(r)
   local idx = Lin.IntVec(r:ndim())

   local ix = r:lower(1)
   for loc = 1, r:volume() do
      invIndexer(loc, idx)
      assert_equal(ix, idx[1], "test_25: Checking inverse indexer")
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
      assert_equal(idx[1], invIdx[1], "test_26: Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "test_26: Checking inv indexer")
   end
end

function test_27()
   local r = Range.Range({-1, -5, 0}, {1, 10, 4})
   local indexer = Range.makeColMajorGenIndexer(r)
   local invIndexer = Range.makeColMajorInvIndexer(r)

   local invIdx = Lin.IntVec(r:ndim())
   for idx in r:colMajorIter() do
      invIndexer(indexer(idx), invIdx)
      assert_equal(idx[1], invIdx[1], "test_27: Checking inv indexer")
      assert_equal(idx[2], invIdx[2], "test_27: Checking inv indexer")
      assert_equal(idx[3], invIdx[3], "test_27: Checking inv indexer")
   end
end

function test_28()
   local r = Range.Range({1, 1}, {10, 10})

   local lowerSkin = r:lowerSkin(1, 1)
   assert_equal(10, lowerSkin:volume(), "test_28: Checking volume")
   assert_equal(1, lowerSkin:lower(1), "test_28: Checking lower in 1")
   assert_equal(1, lowerSkin:upper(1), "test_28: Checking upper in 1")
   assert_equal(1, lowerSkin:lower(2), "test_28: Checking lower in 2")
   assert_equal(10, lowerSkin:upper(2), "test_28: Checking upper in 2")

   local lowerSkin = r:lowerSkin(2, 1)
   assert_equal(10, lowerSkin:volume(), "test_28: Checking volume")
   assert_equal(1, lowerSkin:lower(1), "test_28: Checking lower in 1")
   assert_equal(10, lowerSkin:upper(1), "test_28: Checking upper in 1")
   assert_equal(1, lowerSkin:lower(2), "test_28: Checking lower in 2")
   assert_equal(1, lowerSkin:upper(2), "test_28: Checking upper in 2")

   local upperSkin = r:upperSkin(1, 1)
   assert_equal(10, upperSkin:volume(), "test_28: Checking volume")
   assert_equal(10, upperSkin:lower(1), "test_28: Checking upper in 1")
   assert_equal(10, upperSkin:upper(1), "test_28: Checking upper in 1")
   assert_equal(1, upperSkin:lower(2), "test_28: Checking upper in 2")
   assert_equal(10, upperSkin:upper(2), "test_28: Checking upper in 2")

   local upperSkin = r:upperSkin(2, 1)
   assert_equal(10, upperSkin:volume(), "test_28: Checking volume")
   assert_equal(1, upperSkin:lower(1), "test_28: Checking upper in 1")
   assert_equal(10, upperSkin:upper(1), "test_28: Checking upper in 1")
   assert_equal(10, upperSkin:lower(2), "test_28: Checking upper in 2")
   assert_equal(10, upperSkin:upper(2), "test_28: Checking upper in 2")

   ---

   local lowerGhost = r:lowerGhost(1, 1)
   assert_equal(10, lowerGhost:volume(), "test_28: Checking volume")
   assert_equal(0, lowerGhost:lower(1), "test_28: Checking lower in 1")
   assert_equal(0, lowerGhost:upper(1), "test_28: Checking upper in 1")
   assert_equal(1, lowerGhost:lower(2), "test_28: Checking lower in 2")
   assert_equal(10, lowerGhost:upper(2), "test_28: Checking upper in 2")

   local lowerGhost = r:lowerGhost(2, 1)
   assert_equal(10, lowerGhost:volume(), "test_28: Checking volume")
   assert_equal(1, lowerGhost:lower(1), "test_28: Checking lower in 1")
   assert_equal(10, lowerGhost:upper(1), "test_28: Checking upper in 1")
   assert_equal(0, lowerGhost:lower(2), "test_28: Checking lower in 2")
   assert_equal(0, lowerGhost:upper(2), "test_28: Checking upper in 2")

   local upperGhost = r:upperGhost(1, 1)
   assert_equal(10, upperGhost:volume(), "test_28: Checking volume")
   assert_equal(11, upperGhost:lower(1), "test_28: Checking upper in 1")
   assert_equal(11, upperGhost:upper(1), "test_28: Checking upper in 1")
   assert_equal(1, upperGhost:lower(2), "test_28: Checking upper in 2")
   assert_equal(10, upperGhost:upper(2), "test_28: Checking upper in 2")

   local upperGhost = r:upperGhost(2, 1)
   assert_equal(10, upperGhost:volume(), "test_28: Checking volume")
   assert_equal(1, upperGhost:lower(1), "test_28: Checking upper in 1")
   assert_equal(10, upperGhost:upper(1), "test_28: Checking upper in 1")
   assert_equal(11, upperGhost:lower(2), "test_28: Checking upper in 2")
   assert_equal(11, upperGhost:upper(2), "test_28: Checking upper in 2")

end

function test_29()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})

   local extRange1 = range:extendDirs({1}, 1, 2)
   assert_equal(0, extRange1:lower(1), "test_29: Checking lower extended")
   assert_equal(1, extRange1:lower(2), "test_29: Checking lower extended")
   assert_equal(1, extRange1:lower(3), "test_29: Checking lower extended")

   assert_equal(12, extRange1:upper(1), "test_29: Checking upper extended")
   assert_equal(20, extRange1:upper(2), "test_29: Checking upper extended")
   assert_equal(30, extRange1:upper(3), "test_29: Checking upper extended")

   extRange1 = range:extendDirs({2}, 1, 2)
   assert_equal(1, extRange1:lower(1), "test_29: Checking lower extended")
   assert_equal(0, extRange1:lower(2), "test_29: Checking lower extended")
   assert_equal(1, extRange1:lower(3), "test_29: Checking lower extended")

   assert_equal(10, extRange1:upper(1), "test_29: Checking upper extended")
   assert_equal(22, extRange1:upper(2), "test_29: Checking upper extended")
   assert_equal(30, extRange1:upper(3), "test_29: Checking upper extended")

   extRange1 = range:extendDirs({3}, 1, 2)
   assert_equal(1, extRange1:lower(1), "test_29: Checking lower extended")
   assert_equal(1, extRange1:lower(2), "test_29: Checking lower extended")
   assert_equal(0, extRange1:lower(3), "test_29: Checking lower extended")

   assert_equal(10, extRange1:upper(1), "test_29: Checking upper extended")
   assert_equal(20, extRange1:upper(2), "test_29: Checking upper extended")
   assert_equal(32, extRange1:upper(3), "test_29: Checking upper extended")

   extRange1 = range:extendDirs({1,3}, 1, 2)
   assert_equal(0, extRange1:lower(1), "test_29: Checking lower extended")
   assert_equal(1, extRange1:lower(2), "test_29: Checking lower extended")
   assert_equal(0, extRange1:lower(3), "test_29: Checking lower extended")

   assert_equal(12, extRange1:upper(1), "test_29: Checking upper extended")
   assert_equal(20, extRange1:upper(2), "test_29: Checking upper extended")
   assert_equal(32, extRange1:upper(3), "test_29: Checking upper extended")   
end

function test_30()
   local range = Range.Range({1,1,1,1,1}, {10,20,30,40,50})

   local r1 = range:selectFirst(2)

   assert_equal(2, r1:ndim(), "test_30: Checking dimensions")
   for i = 1, r1:ndim() do
      assert_equal(range:lower(i), r1:lower(i), "test_30: Checking reduced range")
      assert_equal(range:upper(i), r1:upper(i), "test_30: Checking reduced range")
   end

   local r2 = range:selectLast(3)
   assert_equal(3, r2:ndim(), "test_30: Checking dimensions")
   for i = 1, r2:ndim() do
      assert_equal(range:lower(i+2), r2:lower(i), "test_30: Checking reduced range")
      assert_equal(range:upper(i+2), r2:upper(i), "test_30: Checking reduced range")
   end
end

function test_31()
   local range = Range.Range({1,1,1}, {10,20,30})

   -- row-major test
   local idxr = Range.makeGenIndexer(Range.rowMajor, range)
   local lastIdx = 0
   for idx in range:iter(Range.rowMajor) do
      local currIdx = idxr(idx)
      assert_equal(1, currIdx-lastIdx, "test_31: Checking generic indexer/iterator pair")
      lastIdx = currIdx
   end

   -- col-major test
   idxr = Range.makeGenIndexer(Range.colMajor, range)
   lastIdx = 0
   for idx in range:iter(Range.colMajor) do
      local currIdx = idxr(idx)
      assert_equal(1, currIdx-lastIdx, "test_31: Checking generic indexer/iterator pair")
      lastIdx = currIdx
   end   
end

function test_32()
   local range = Range.Range({1,1,1}, {10,20,30})

   -- row-major test
   local idxr = range:genIndexer(Range.rowMajor)
   local lastIdx = 0
   for idx in range:iter(Range.rowMajor) do
      local currIdx = idxr(idx)
      assert_equal(1, currIdx-lastIdx, "test_32: Checking generic indexer/iterator pair")
      lastIdx = currIdx
   end

   -- col-major test
   idxr = range:genIndexer(Range.colMajor)
   lastIdx = 0
   for idx in range:iter(Range.colMajor) do
      local currIdx = idxr(idx)
      assert_equal(1, currIdx-lastIdx, "test_32: Checking generic indexer/iterator pair")
      lastIdx = currIdx
   end   
end

function test_33()
   local range = Range.Range({1,1}, {10,20})

   -- row-major test
   local idxr = range:indexer(Range.rowMajor)
   local lastIdx = 0

   for i = range:lower(1), range:upper(1) do
      for j = range:lower(2), range:upper(2) do
	 local currIdx = idxr(i,j)
	 assert_equal(1, currIdx-lastIdx, "test_33: Checking generic indexer/iterator pair")
	 lastIdx = currIdx
      end
   end

   -- col-major test
   idxr = range:indexer(Range.colMajor)
   lastIdx = 0
   for j = range:lower(2), range:upper(2) do   
      for i = range:lower(1), range:upper(1) do
	 local currIdx = idxr(i, j)
	 assert_equal(1, currIdx-lastIdx, "test_33: Checking generic indexer/iterator pair")
	 lastIdx = currIdx
      end
   end   
end

function test_34()
   if not GKYL_HAVE_CUDA then return end

   local range = Range.Range({0, 0}, {1, 5})
   local cuRange, err = Range.copyHostToDevice(range)
   assert_equal(0, err, "test_34: Checking if range object copied to device")
end

function test_35()
   local range = Range.Range({1,1,1}, {10,20,30})
   assert_equal(true, range:contains({1,1,1}), "test_35: Checking contains")
   assert_equal(true, range:contains({10,20,30}), "test_35: Checking contains")
   assert_equal(true, range:contains({3,3,3}), "test_35: Checking contains")   

   assert_equal(false, range:contains({0,1,1}), "test_35: Checking contains")
   assert_equal(false, range:contains({10,21,30}), "test_35: Checking contains")
		
end

function test_36()
   if not GKYL_HAVE_CUDA then return end
   local range = Range.Range({0, 0}, {1, 5})
   local devRange = range:cloneOnDevice()
   ffi.C.unit_showRange(devRange)   
end

function test_37()
   local r1 = Range.Range({1, 2}, {10, 11})

   local r2 = Range.Range({5, 6}, {12, 13})
   local r3 = r1:difference(r2)
   assert_equal(1, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(2, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(4, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(5, r3:upper(2), "test_37: Checking range relative difference")

   local r3 = r2:difference(r1)
   assert_equal(11, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(12, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(12, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(13, r3:upper(2), "test_37: Checking range relative difference")

   local r2 = Range.Range({5, 5}, {10, 10})
   local r3 = r1:difference(r2)
   assert_equal(1, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(2, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(4, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(4, r3:upper(2), "test_37: Checking range relative difference")

   local r2 = Range.Range({4, 4}, {15, 8})
   local r3 = r1:difference(r2)
   assert_equal(1, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(2, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(3, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(3, r3:upper(2), "test_37: Checking range relative difference")

   local r2 = Range.Range({1, 2}, {10, 11})
   assert_equal(true, r1:isDifferenceEmpty(r2), "test_37: Checking range relative difference")

   local r2 = Range.Range({5, 2}, {16, 11})
   local r3 = r1:difference(r2)
   assert_equal(1, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(2, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(4, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(11, r3:upper(2), "test_37: Checking range relative difference")

   local r2 = Range.Range({14, 15}, {21, 22})
   local r3 = r1:difference(r2)
   assert_equal(1, r3:lower(1), "test_37: Checking range relative difference")
   assert_equal(2, r3:lower(2), "test_37: Checking range relative difference")
   assert_equal(10, r3:upper(1), "test_37: Checking range relative difference")
   assert_equal(11, r3:upper(2), "test_37: Checking range relative difference")
   assert_equal(false, r1:isDifferenceEmpty(r2), "test_37: Checking range relative difference")

end

function test_38()
   local r = Range.Range({1, 2}, {10, 11})
   local sub = r:subRange({3, 5}, {7, 9})

   assert_equal(1, sub:isSubRange(), "test_38: Checking isSubRange")
end

function test_39()
   local r = Range.Range({1, 2}, {10, 11})
   local ext = r:extend(1, 1)
   local sub = ext:subRange(r:lowerAsVec(), r:upperAsVec())

   assert_equal(1, sub:isSubRange(), "test_39: Checking isSubRange")
   assert_equal(0, ext:isSubRange(), "test_39: Checking isSubRange")
   assert_equal(0, r:isSubRange(), "test_39: Checking isSubRange")
   
   assert_equal(1, r.nsplit)
   assert_equal(1, ext.nsplit)
   assert_equal(1, sub.nsplit)

   assert_equal(r:lower(1), sub:lower(1))
   assert_equal(r:lower(2), sub:lower(2))
   assert_equal(r:upper(1), sub:upper(1))
   assert_equal(r:upper(2), sub:upper(2))
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
test_15_b()
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
test_30()
test_31()
test_32()
test_33()
test_34()
test_35()
if GKYL_HAVE_CUDA then
   test_36()
end
test_37()
test_38()
test_39()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
