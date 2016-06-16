-- Gkyl ------------------------------------------------------------------------
--
-- Test for fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
   }

   assert_equal(1, field:ndim(), "Checking dimensions")
   assert_equal(1, field:lowerGhost(), "Checking lower ghost")
   assert_equal(1, field:upperGhost(), "Checking upper ghost")
   assert_equal((10+2)*3, field:size(), "Checking size")

   assert_equal("col-major", field:layout(), "Checking layout")

   local localRange = field:localRange()
   assert_equal(1, localRange:lower(1), "Checking range lower")
   assert_equal(10, localRange:upper(1), "Checking range upper")

   local localExtRange = field:localExtRange()
   assert_equal(0, localExtRange:lower(1), "Checking range lower")
   assert_equal(11, localExtRange:upper(1), "Checking range upper")

   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end

   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      assert_equal(i+1, fitr[1], "Checking field value")
      assert_equal(i+2, fitr[2], "Checking field value")
      assert_equal(i+3, fitr[3], "Checking field value")
   end
end

function test_2()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }

   assert_equal(2, field:ndim(), "Checking dimensions")
   assert_equal(1, field:lowerGhost(), "Checking lower ghost")
   assert_equal(2, field:upperGhost(), "Checking upper ghost")
   assert_equal((10+3)*(10+3)*3, field:size(), "Checking size")

   assert_equal("col-major", field:layout(), "Checking layout")

   local localRange = field:localRange()
   assert_equal(1, localRange:lower(1), "Checking range lower")
   assert_equal(1, localRange:lower(2), "Checking range lower")

   assert_equal(10, localRange:upper(1), "Checking range upper")
   assert_equal(10, localRange:upper(2), "Checking range upper")

   local localExtRange = field:localExtRange()
   assert_equal(0, localExtRange:lower(1), "Checking range lower")
   assert_equal(0, localExtRange:lower(2), "Checking range lower")

   assert_equal(12, localExtRange:upper(1), "Checking range upper")
   assert_equal(12, localExtRange:upper(2), "Checking range upper")

   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr = field:get(indexer(i,j))
	 fitr[1] = i+2*j+1
	 fitr[2] = i+2*j+2
	 fitr[3] = i+2*j+3
      end
   end

   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr = field:get(indexer(i,j))
	 assert_equal(i+2*j+1, fitr[1], "Checking field value")
	 assert_equal(i+2*j+2, fitr[2], "Checking field value")
	 assert_equal(i+2*j+3, fitr[3], "Checking field value")
      end
   end
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local EulerField = DataStruct.new_field_ct(ffi.typeof("struct {double rho, rhou, E;}"))
   local field = EulerField {
      onGrid = grid,
      ghost = {1, 1},
   }

   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1].rho = 1
      fitr[1].rhou = 2
      fitr[1].E = 3
   end   
end

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
