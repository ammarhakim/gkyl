-- Gkyl ------------------------------------------------------------------------
--
-- Test for fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi        = require "ffi"
local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
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
   for i = localExtRange:lower(1), localExtRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end

   for j = 1, field:numComponents() do
      for i = localExtRange:lower(1), localExtRange:upper(1) do
         local idx = indexer(i)
         assert_equal(i+j, field._data[(j-1)+3*(idx-1)], "Checking values by indexing _data directly")
         assert_equal(i+j, field:get(idx)[j], "Checking values by using get()")         
         assert_equal(i+j, field:getDataPtrAt(idx)[j-1], "Checking values by using getDataPtrAt()")         
      end
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
      createDeviceCopy = false,
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
      metaData = {
	 polyOrder = 2, basisType = "ms"
      }
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

function test_13()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      periodicDirs = {1, 2},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }
   field:clear(10.5)

   -- Set corner cells.
   local indexer = field:indexer()
   local fItr
   fItr = field:get(indexer(1,1));   fItr[1] = 1.0
   fItr = field:get(indexer(10,1));  fItr[1] = 2.0
   fItr = field:get(indexer(1,10));  fItr[1] = 3.0
   fItr = field:get(indexer(10,10)); fItr[1] = 4.0

   field:sync() -- Sync field.

   -- Check if periodic dirs are sync()-ed properly.
   local fItr = field:get(indexer(11,1))
   assert_equal(1.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(11,10))
   assert_equal(3.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,1))
   assert_equal(2.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,10))
   assert_equal(4.0, fItr[1], "Checking non-corner periodic sync")

   -- Check corner ghost cells.
   local fItr = field:get(indexer(11,11))
   assert_equal(1.0, fItr[1], "Checking 11,11 corner periodic sync")
   local fItr = field:get(indexer(11,0))
   assert_equal(3.0, fItr[1], "Checking 11,0 corner periodic sync")
   local fItr = field:get(indexer(0,0))
   assert_equal(4.0, fItr[1], "Checking 0,0 corner periodic sync")
   local fItr = field:get(indexer(0,11))
   assert_equal(2.0, fItr[1], "Checking 0,11 corner periodic sync")
end

function test_14()
   -- Test syncCorners in 3D.
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {1.0, 1.0, 1.0},
      cells = {10, 10, 10},
      periodicDirs = {1, 2, 3},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }
   field:clear(10.5)

   local indexer = field:indexer()
   local fItr
   -- Set corner cells.
   fItr = field:get(indexer(1,1,1));   fItr[1] = 1.0
   fItr = field:get(indexer(10,1,1));  fItr[1] = 2.0
   fItr = field:get(indexer(1,10,1));  fItr[1] = 3.0
   fItr = field:get(indexer(10,10,1)); fItr[1] = 4.0
   fItr = field:get(indexer(1,1,8));   fItr[1] = 5.0
   fItr = field:get(indexer(10,1,8));  fItr[1] = 6.0
   fItr = field:get(indexer(1,10,8));  fItr[1] = 7.0
   fItr = field:get(indexer(10,10,8)); fItr[1] = 8.0

   -- Set skin cells that aren't corners. 
   for i = 2,9 do fItr = field:get(indexer(i,1,1));   fItr[1] = 1.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(1,i,1));   fItr[1] = 2.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(i,10,1));  fItr[1] = 3.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(10,i,1));  fItr[1] = 4.0+i/10. end

   for i = 2,9 do fItr = field:get(indexer(i,1,10));  fItr[1] = 5.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(1,i,10));  fItr[1] = 6.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(i,10,10)); fItr[1] = 7.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(10,i,10)); fItr[1] = 8.0+i/10. end

   for i = 2,9 do fItr = field:get(indexer(1,1,i));   fItr[1] = 11.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(10,1,i));  fItr[1] = 12.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(1,10,i));  fItr[1] = 13.0+i/10. end
   for i = 2,9 do fItr = field:get(indexer(10,10,i)); fItr[1] = 14.0+i/10. end

   field:sync() -- Sync field.

   -- Check if periodic dirs are sync()-ed properly.
   local fItr = field:get(indexer(11,1,1))
   assert_equal(1.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(11,10,1))
   assert_equal(3.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,1,1))
   assert_equal(2.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,10,1))
   assert_equal(4.0, fItr[1], "Checking non-corner periodic sync")

   -- Check corner ghost cells.
   local fItr = field:get(indexer(11,11,11))
   assert_equal(1.0, fItr[1], "Checking 11,11,11 corner periodic sync")
   local fItr = field:get(indexer(11,0,11))
   assert_equal(3.0, fItr[1], "Checking 11,0,11 corner periodic sync")
   local fItr = field:get(indexer(0,0,11))
   assert_equal(4.0, fItr[1], "Checking 0,0,11 corner periodic sync")
   local fItr = field:get(indexer(0,11,11))
   assert_equal(2.0, fItr[1], "Checking 0,11,11 corner periodic sync")
   -- Check "corner" ranges (edges really), which are not actual corners of the extended range.
   for i = 2,9 do
      fItr = field:get(indexer(i,11,11))
      assert_equal(1.0+i/10., fItr[1], string.format("Checking %d,11,11 corner periodic sync",i))
      fItr = field:get(indexer(11,i,11))
      assert_equal(2.0+i/10., fItr[1], string.format("Checking 11,%d,11 corner periodic sync",i))
      fItr = field:get(indexer(i,0,11))
      assert_equal(3.0+i/10., fItr[1], string.format("Checking %d,0,11 corner periodic sync",i))
      fItr = field:get(indexer(0,i,11))
      assert_equal(4.0+i/10., fItr[1], string.format("Checking 0,%d,11 corner periodic sync",i))

      fItr = field:get(indexer(i,11,0))
      assert_equal(5.0+i/10., fItr[1], string.format("Checking %d,11,0 corner periodic sync",i))
      fItr = field:get(indexer(11,i,0))
      assert_equal(6.0+i/10., fItr[1], string.format("Checking 11,%d,0 corner periodic sync",i))
      fItr = field:get(indexer(i,0,0))
      assert_equal(7.0+i/10., fItr[1], string.format("Checking %d,0,0 corner periodic sync",i))
      fItr = field:get(indexer(0,i,0))
      assert_equal(8.0+i/10., fItr[1], string.format("Checking 0,%d,0 corner periodic sync",i))

      fItr = field:get(indexer(11,11,i))
      assert_equal(11.0+i/10., fItr[1], string.format("Checking 11,11,%d corner periodic sync",i))
      fItr = field:get(indexer(0,11,i))
      assert_equal(12.0+i/10., fItr[1], string.format("Checking 0,11,%d corner periodic sync",i))
      fItr = field:get(indexer(11,0,i))
      assert_equal(13.0+i/10., fItr[1], string.format("Checking 11,0,%d corner periodic sync",i))
      fItr = field:get(indexer(0,0,i))
      assert_equal(14.0+i/10., fItr[1], string.format("Checking 0,0,%d corner periodic sync",i))
   end
end

function test_15()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      periodicDirs = {1, 2},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }
   field:clear(10.5)

   -- Set corner cells.
   local indexer = field:indexer()
   local fItr
   fItr = field:get(indexer(1,1));   fItr[1] = 1.0
   fItr = field:get(indexer(10,1));  fItr[1] = 2.0
   fItr = field:get(indexer(1,10));  fItr[1] = 3.0
   fItr = field:get(indexer(10,10)); fItr[1] = 4.0

   field:periodicCopy() -- Copy periodic boundary conditions for field. SyncCorners not implemented for this method yet.

   -- Check if periodic dirs are sync()-ed properly.
   local fItr = field:get(indexer(11,1))
   assert_equal(1.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(11,10))
   assert_equal(3.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,1))
   assert_equal(2.0, fItr[1], "Checking non-corner periodic sync")
   local fItr = field:get(indexer(0,10))
   assert_equal(4.0, fItr[1], "Checking non-corner periodic sync")

end

function test_16()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 2},
   }
   local scalar = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 2},
   }
   field:clear(10.0)
   scalar:clear(2.5)
   field:scaleByCell(scalar)

   local indexer = field:genIndexer()
   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(25.0, fitr[1], "Checking scaled field value")
      assert_equal(25.0, fitr[2], "Checking scaled field value")
      assert_equal(25.0, fitr[3], "Checking scaled field value")
   end   

   for idx in scalar:localExtRangeIter() do
      local sitr = scalar:get(indexer(idx))
      sitr[1] = idx[1]
   end   

   field:scaleByCell(scalar)

   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(25.0*idx[1], fitr[1], "Checking scaled field value")
      assert_equal(25.0*idx[1], fitr[2], "Checking scaled field value")
      assert_equal(25.0*idx[1], fitr[3], "Checking scaled field value")
   end   

   -- Initialize fields to random numbers.
   math.randomseed(1000*os.time())
   local localRange = scalar:localRange()
   local fldIdxr    = field:genIndexer()
   local scaIdxr    = scalar:genIndexer()
   for idx in localRange:rowMajorIter() do
      local fldItr = field:get(fldIdxr( idx ))
      local scaItr = scalar:get(scaIdxr( idx ))
      for k = 1, field:numComponents() do
         fldItr[k] = math.random()
      end
      for k = 1, scalar:numComponents() do
         scaItr[k] = math.random()
      end
   end
   -- Get the maximum by stepping through the scalar.
   fldMax, fldMin, fldSum =  {}, {}, {}
   for k = 1, field:numComponents() do
      fldMax[k] = GKYL_MIN_DOUBLE
      fldMin[k] = GKYL_MAX_DOUBLE
      fldSum[k] = 0.0
   end
   scaMax, scaMin, scaSum = GKYL_MIN_DOUBLE, GKYL_MAX_DOUBLE, 0.0
   for idx in localRange:rowMajorIter() do
      local fldItr = field:get(fldIdxr( idx ))
      local scaItr = scalar:get(scaIdxr( idx ))
      for k = 1, field:numComponents() do
         fldMax[k] = math.max(fldMax[k],fldItr[k])
         fldMin[k] = math.min(fldMin[k],fldItr[k])
         fldSum[k] = fldSum[k] + fldItr[k]
      end
      scaMax = math.max(scaMax,scaItr[1])
      scaMin = math.min(scaMin,scaItr[1])
      scaSum = scaSum + scaItr[1]
   end
   cartFldMax, cartFldMin, cartFldSum = field:reduce("max"), field:reduce("min"), field:reduce("sum")
   cartScaMax, cartScaMin, cartScaSum = scalar:reduce("max"), scalar:reduce("min"), scalar:reduce("sum")
   for k = 1, field:numComponents() do
      assert_equal(fldMax[k], cartFldMax[k], "Checking reduce('max')")
      assert_equal(fldMin[k], cartFldMin[k], "Checking reduce('min')")
      assert_equal(fldSum[k], cartFldSum[k], "Checking reduce('sum')")
   end
   assert_equal(scaMax, cartScaMax[1], "Checking scalar reduce('max')")
   assert_equal(scaMin, cartScaMin[1], "Checking scalar reduce('min')")
   assert_equal(scaSum, cartScaSum[1], "Checking scalar reduce('sum')")
end

function test_17()
   -- Test the :accumulateOffset and :combineOffset methods.
   -- NOTE: the offset in these methods are field component offsets, not vector offsets.
   --       For a CartField with multiple DG fields, selecting the ith DG field would
   --       require the offset (i-1)*numBasis.
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field6 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 6,
      ghost         = {1, 2},
   }
   local field8 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 8,
      ghost         = {1, 2},
   }
   field6:clear(20.0)
   field8:clear(10.0)

   -- field1 and field3 have fewer components than field8.
   local field1 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 2},
   }
   local field3 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 2},
   }
   field1:clear(0.0)
   field3:clear(0.0)

   local indexer = field1:genIndexer()
   for idx in field1:localExtRangeIter() do
      local fitr = field1:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
   end

   local indexer = field3:genIndexer()
   for idx in field3:localExtRangeIter() do
      local fitr = field3:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   -- Accumulate and combine onto field with more components.
   field8:accumulateOffset(1.0, field3, 0, 2.0, field3, 5)
   field6:combineOffset(1.0, field1, 0, 2.0, field1, 3)

   local indexer = field8:genIndexer()
   for idx in field8:localExtRangeIter() do
      local fitr = field8:get(indexer(idx))
      assert_equal(10+(idx[1]+2*idx[2]+1), fitr[1], "Checking field8 value")
      assert_equal(10+(idx[1]+2*idx[2]+2), fitr[2], "Checking field8 value")
      assert_equal(10+(idx[1]+2*idx[2]+3), fitr[3], "Checking field8 value")
      assert_equal(10, fitr[4], "Checking field8 value")
      assert_equal(10, fitr[5], "Checking field8 value")
      assert_equal(10+2*(idx[1]+2*idx[2]+1), fitr[6], "Checking field8 value")
      assert_equal(10+2*(idx[1]+2*idx[2]+2), fitr[7], "Checking field8 value")
      assert_equal(10+2*(idx[1]+2*idx[2]+3), fitr[8], "Checking field8 value")
   end

   local indexer = field6:genIndexer()
   for idx in field6:localExtRangeIter() do
      local fitr = field6:get(indexer(idx))
      assert_equal((idx[1]+2*idx[2]+1), fitr[1], "Checking field6 value")
      assert_equal(20, fitr[2], "Checking field6 value")
      assert_equal(20, fitr[3], "Checking field6 value")
      assert_equal(2*(idx[1]+2*idx[2]+1), fitr[4], "Checking field6 value")
      assert_equal(20, fitr[5], "Checking field6 value")
      assert_equal(20, fitr[6], "Checking field6 value")
   end

   -- Test accumulating and combining onto field with fewer components.
   field1:clear(20.0)
   field3:clear(10.0)
   field6:clear(0.0)
   field8:clear(0.0)

   local indexer = field8:genIndexer()
   for idx in field8:localExtRangeIter() do
      local fitr = field8:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
      fitr[4] = idx[1]+2*idx[2]+4
      fitr[5] = idx[1]+2*idx[2]+5
      fitr[6] = idx[1]+2*idx[2]+6
      fitr[7] = idx[1]+2*idx[2]+7
      fitr[8] = idx[1]+2*idx[2]+8
   end

   local indexer = field6:genIndexer()
   for idx in field6:localExtRangeIter() do
      local fitr = field6:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
      fitr[4] = idx[1]+2*idx[2]+4
      fitr[5] = idx[1]+2*idx[2]+5
      fitr[6] = idx[1]+2*idx[2]+6
   end

   field3:accumulateOffset(1.0, field8, 2, 2.0, field8, 5)
   field1:combineOffset(1.0, field6, 1, 2.0, field6, 4)

   local indexer = field3:genIndexer()
   for idx in field3:localExtRangeIter() do
      local fitr = field3:get(indexer(idx))
      assert_equal(10+(idx[1]+2*idx[2]+3)+2.0*(idx[1]+2*idx[2]+6), fitr[1], "Checking field3 value")
      assert_equal(10+(idx[1]+2*idx[2]+4)+2.0*(idx[1]+2*idx[2]+7), fitr[2], "Checking field3 value")
      assert_equal(10+(idx[1]+2*idx[2]+5)+2.0*(idx[1]+2*idx[2]+8), fitr[3], "Checking field3 value")
   end

   local indexer = field1:genIndexer()
   for idx in field1:localExtRangeIter() do
      local fitr = field1:get(indexer(idx))
      assert_equal((idx[1]+2*idx[2]+2)+2.0*(idx[1]+2*idx[2]+5), fitr[1], "Checking field1 value")
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
test_13()
test_14()
test_15()
test_16()
test_17()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
