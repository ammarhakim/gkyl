-- Gkyl ------------------------------------------------------------------------
--
-- Test to illustrate looping using linear index into a box
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Range = require "Lib.Range"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {20, 20},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1},
   }
   field:clear(0.0)
   
   local globalRange = grid:globalRange()
   local invIndexer = Range.makeRowMajorInvIndexer(globalRange)
   local idx = Lin.IntVec(globalRange:ndim())

   -- sub-range to index over
   local localRange = Range.Range ( {3, 5}, {12, 14} )
   local indexer = field:genIndexer()
   
   -- loop over global range but only update local-range
   for loc = 1, globalRange:volume() do
      invIndexer(loc, idx)
      -- check if this index is in the local range
      if localRange:contains(idx) then
	 -- perform update
	 local fitr = field:get(indexer(idx))
	 fitr[1] = 10.5
      end
   end

   -- check if setting worked
   for idx in globalRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      if localRange:contains(idx) then
	 assert_equal(10.5, fitr[1], "Checking interior")
      else
	 assert_equal(0.0, fitr[1], "Checking exterior")
      end
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
