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

   assert_equal(field:defaultLayout(), field:layout(), "Checking layout")

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

   assert_equal(field:defaultLayout(), field:layout(), "Checking layout")

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
      fitr[i].rho = i+1
      fitr[i].rhou = i+2
      fitr[i].E = i+3
   end

   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      assert_equal(i+1, fitr[i].rho, "Testing Field of struct")
      assert_equal(i+2, fitr[i].rhou, "Testing Field of struct")
      assert_equal(i+3, fitr[i].E, "Testing Field of struct")
   end
   
end

function test_4()
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

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], "Checking field value")
      assert_equal(idx[1]+2*idx[2]+2, fitr[2], "Checking field value")
      assert_equal(idx[1]+2*idx[2]+3, fitr[3], "Checking field value")
   end
end

function test_5()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field1 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
   }
   local field2 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
   }

   local localRange, indexer = field1:localRange(), field1:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field1:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end
   -- copy and test
   field2:copy(field1)
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr1 = field1:get(indexer(i))
      local fitr2 = field2:get(indexer(i))

      assert_equal(fitr1[1], fitr2[1], "Checking if copy worked")
      assert_equal(fitr1[2], fitr2[2], "Checking if copy worked")
      assert_equal(fitr1[3], fitr2[3], "Checking if copy worked")
   end   
   
end

function test_6()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field1 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
   }

   local field2 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
   }

   local fitr = field1:get(1)
   local localRange, indexer = field1:localRange(), field1:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      field1:fill(indexer(i), fitr)
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end
   -- copy and test
   field2:copy(field1)
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr1 = field1:get(indexer(i))
      local fitr2 = field2:get(indexer(i))

      assert_equal(fitr1[1], fitr2[1], "Checking if copy worked")
      assert_equal(fitr1[2], fitr2[2], "Checking if copy worked")
      assert_equal(fitr1[3], fitr2[3], "Checking if copy worked")
   end
end

function test_7()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
   }

   local localRange, indexer = field:localExtRange(), field:indexer()

   -- clear it
   field:clear(1.5)
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))

      assert_equal(1.5, fitr[1], "Checking if clear worked")
      assert_equal(1.5, fitr[2], "Checking if clear worked")
      assert_equal(1.5, fitr[3], "Checking if clear worked")
   end

   field:clear(2.5)
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))

      assert_equal(2.5, fitr[1], "Checking if clear worked")
      assert_equal(2.5, fitr[2], "Checking if clear worked")
      assert_equal(2.5, fitr[3], "Checking if clear worked")
   end

   field:clear(0.0)
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))

      assert_equal(0.0, fitr[1], "Checking if clear worked")
      assert_equal(0.0, fitr[2], "Checking if clear worked")
      assert_equal(0.0, fitr[3], "Checking if clear worked")
   end
end

function test_8()
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

   local indexer = field:genIndexer()
   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], "Checking field value")
      assert_equal(idx[1]+2*idx[2]+2, fitr[2], "Checking field value")
      assert_equal(idx[1]+2*idx[2]+3, fitr[3], "Checking field value")
   end
end

function test_9()
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
   field:clear(10.0)

   local field1 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }

   local indexer = field1:genIndexer()
   for idx in field1:localExtRangeIter() do
      local fitr = field1:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   -- accumulate stuff
   field:accumulate(1.0, field1, 2.0, field1)

   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(10+3*(idx[1]+2*idx[2]+1), fitr[1], "Checking field value")
      assert_equal(10+3*(idx[1]+2*idx[2]+2), fitr[2], "Checking field value")
      assert_equal(10+3*(idx[1]+2*idx[2]+3), fitr[3], "Checking field value")
   end

   -- combine stuff
   field:combine(1.0, field1, 2.0, field1)

   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(3*(idx[1]+2*idx[2]+1), fitr[1], "Checking field value")
      assert_equal(3*(idx[1]+2*idx[2]+2), fitr[2], "Checking field value")
      assert_equal(3*(idx[1]+2*idx[2]+3), fitr[3], "Checking field value")
   end   
end

function test_10()
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
   field:clear(10.0)
   field:scale(2.5)

   local indexer = field:genIndexer()
   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(25.0, fitr[1], "Checking scaled field value")
      assert_equal(25.0, fitr[2], "Checking scaled field value")
      assert_equal(25.0, fitr[3], "Checking scaled field value")
   end   
end

function test_11()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3*3,
      ghost = {1, 2},
   }

   field:clear(10.0)

   local function set(fIn, val)
      for i = 1, 3 do fIn[i] = val end
   end

   local indexer = field:genIndexer()
   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      local data = fitr:data()
      set(data+0-1, 1.5)
      set(data+3-1, 2.5)
      set(data+6-1, 3.5)
   end

   local indexer = field:genIndexer()
   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(1.5, fitr[1], "Checking set")
      assert_equal(1.5, fitr[2], "Checking set")
      assert_equal(1.5, fitr[3], "Checking set")

      assert_equal(2.5, fitr[4], "Checking set")
      assert_equal(2.5, fitr[5], "Checking set")
      assert_equal(2.5, fitr[6], "Checking set")

      assert_equal(3.5, fitr[7], "Checking set")
      assert_equal(3.5, fitr[8], "Checking set")
      assert_equal(3.5, fitr[9], "Checking set")
   end   
   
end

function test_12()
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
   field:clear(10.25)

   -- write field
   field:write("CartFieldTest_field.bp", 2.5, 50)

   local fieldIn = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }
   fieldIn:clear(0.0)

   local tm, fr = fieldIn:read("CartFieldTest_field.bp")

   assert_equal(2.5, tm, "Checking time-stamp")
   assert_equal(50, fr, "Checking frame")
   
   -- check if fields are identical
   local indexer = field:genIndexer()
   for idx in field:localRangeIter() do
      local fitr, fitrIn = field:get(indexer(idx)), fieldIn:get(indexer(idx))
      for k = 1, field:numComponents() do
	 assert_equal(fitr[k], fitrIn[k], "Checking if field read correctly")
      end
   end
end

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

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
