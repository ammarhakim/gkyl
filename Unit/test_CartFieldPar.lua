-- Gkyl ------------------------------------------------------------------------
--
-- Test for fields on cartesian grids (parallel).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi              = require "ffi"
local Unit             = require "Unit"
local Grid             = require "Grid"
local DecompRegionCalc = require "Lib.CartDecomp"
local DataStruct       = require "DataStruct"
local Mpi              = require "Comm.Mpi"

local cuda
local cuAlloc
if GKYL_HAVE_CUDA then
   cuda    = require "Cuda.RunTime"
   cuAlloc = require "Cuda.Alloc"
end

local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then
      print(msg)
   end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

function test_1(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 2 then
      log("Not running test_1 as numProcs not exactly 2")
      return
   end
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2} }
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }

   assert_equal(1, field:ndim(), "Checking dimensions")
   assert_equal(1, field:lowerGhost(), "Checking lower ghost")
   assert_equal(1, field:upperGhost(), "Checking upper ghost")
   assert_equal(field:localExtRange():volume()*3, field:size(), "Checking size")

   assert_equal(field:defaultLayout(), field:layout(), "Checking layout")

   local localRange = field:localRange()
   if rank == 0 then
      assert_equal(1, localRange:lower(1), "Checking range lower")
      assert_equal(5, localRange:upper(1), "Checking range upper")
   elseif rank == 1 then
      assert_equal(6, localRange:lower(1), "Checking range lower")
      assert_equal(10, localRange:upper(1), "Checking range upper")
   end

   local localExtRange = field:localExtRange()
   if rank == 0 then
      assert_equal(0, localExtRange:lower(1), "Checking ext range lower")
      assert_equal(6, localExtRange:upper(1), "Checking ext range upper")
   elseif rank == 1 then
      assert_equal(5, localExtRange:lower(1), "Checking ext range lower")
      assert_equal(11, localExtRange:upper(1), "Checking ext range upper")
   end

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
   Mpi.Barrier(comm)
end

function test_2(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 4 then
      log("Not running test_2 as numProcs not exactly 4")
      return
   end
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }

   assert_equal(2, field:ndim(), "Checking dimensions")
   assert_equal(1, field:lowerGhost(), "Checking lower ghost")
   assert_equal(2, field:upperGhost(), "Checking upper ghost")
   assert_equal(field:localExtRange():volume()*3, field:size(), "Checking size")

   assert_equal(field:defaultLayout(), field:layout(), "Checking layout")

   local localRange = field:localRange()
   local globalRange = field:globalRange()
   local localExtRange = field:localExtRange()

   assert_equal(100, globalRange:volume(), "Testing global volume")
   assert_equal(25, localRange:volume(), "Testing local volume")
   assert_equal(64, localExtRange:volume(), "Testing local ext volume")

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
   Mpi.Barrier(comm)
end

function test_3(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 3 then
      log("Not running test_3 as numProcs not exactly 3")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {3} }
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {12},
      decomposition = decomp,
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
      fitr[1] = idx[1]+1
      fitr[2] = idx[1]+2
      fitr[3] = idx[1]+3
   end

   field:sync() -- sync up ghost cells

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+1, fitr[1], "Checking field value")
      assert_equal(idx[1]+2, fitr[2], "Checking field value")
      assert_equal(idx[1]+3, fitr[3], "Checking field value")
   end

   Mpi.Barrier(comm)
end

function test_4(comm)
   -- Check sync. Check the corner ghost cells which lie in the center of the domain too.
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 4 then
      log("Not running test_4 as numProcs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }

   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 2},
      syncCorners   = true,
   }

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   field:sync()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
      assert_equal(idx[1]+2*idx[2]+2, fitr[2], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
      assert_equal(idx[1]+2*idx[2]+3, fitr[3], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)
end

function test_5(comm)
   -- Check sync. Even though syncCorners=true here, the corners are not
   -- synced because of the decomposition used below. In this case the corner
   -- ghost cells would be filled by BCs.
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 2 then
      log("Not running test_5 as numProcs not exactly 2")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1} }

   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {500, 500},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }

   local localRange = field:localRange()
   local indexer    = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1]    = idx[1]+2*idx[2]+1
   end

   field:sync()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)
end

function test_6(comm)
   -- Check sync, and the corners shared by the two MPI ranks.
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 2 then
      log("Not running test_5 as numProcs not exactly 2")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2,1} }

   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {6, 8},
      decomposition = decomp,
      periodicDirs  = {1,2},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }

   local localRange = field:localRange()
   local indexer    = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1]    = idx[1]+2*idx[2]+1
   end

   field:sync()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)

   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

   -- Check corner ghost cells (cuts={2,1}, periodicDirs={1,2}).
   if rank==0 then
      local fItr = field:get(indexer({0,0}))
      assert_equal(6+2*8+1, fItr[1], "Checking 0,0 corner periodic sync")
      local fItr = field:get(indexer({4,0}))
      assert_equal(4+2*8+1, fItr[1], "Checking 4,0 corner periodic sync")
      local fItr = field:get(indexer({0,9}))
      assert_equal(6+2*1+1, fItr[1], "Checking 0,9 corner periodic sync")
      local fItr = field:get(indexer({4,9}))
      assert_equal(4+2*1+1, fItr[1], "Checking 4,9 corner periodic sync")
   else
      local fItr = field:get(indexer({3,0}))
      assert_equal(3+2*8+1, fItr[1], "Checking 3,0 corner periodic sync")
      local fItr = field:get(indexer({7,0}))
      assert_equal(1+2*8+1, fItr[1], "Checking 7,0 corner periodic sync")
      local fItr = field:get(indexer({3,9}))
      assert_equal(3+2*1+1, fItr[1], "Checking 3,9 corner periodic sync")
      local fItr = field:get(indexer({7,9}))
      assert_equal(1+2*1+1, fItr[1], "Checking 7,9 corner periodic sync")
   end

--   -- Check corner ghost cells (cuts={1,2}, periodicDirs={1,2}).
--   if rank==0 then
--      local fItr = field:get(indexer({0,0}))
--      assert_equal(6+2*8+1, fItr[1], "Checking 0,0 corner periodic sync")
--      local fItr = field:get(indexer({7,0}))
--      assert_equal(1+2*8+1, fItr[1], "Checking 7,0 corner periodic sync")
--      local fItr = field:get(indexer({0,5}))
--      assert_equal(6+2*5+1, fItr[1], "Checking 0,5 corner periodic sync")
--      local fItr = field:get(indexer({7,5}))
--      assert_equal(1+2*5+1, fItr[1], "Checking 7,5 corner periodic sync")
--   else
--      local fItr = field:get(indexer({0,4}))
--      assert_equal(6+2*4+1, fItr[1], "Checking 0,4 corner periodic sync")
--      local fItr = field:get(indexer({7,4}))
--      assert_equal(1+2*4+1, fItr[1], "Checking 7,4 corner periodic sync")
--      local fItr = field:get(indexer({0,9}))
--      assert_equal(6+2*1+1, fItr[1], "Checking 0,9 corner periodic sync")
--      local fItr = field:get(indexer({7,9}))
--      assert_equal(1+2*1+1, fItr[1], "Checking 7,9 corner periodic sync")
--   end

--   -- Check corner ghost cells (cuts={2,1}, periodicDirs={2}).
--   if rank==0 then
--      local fItr = field:get(indexer({0,0}))
--      assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
--      local fItr = field:get(indexer({4,0}))
--      assert_equal(4+2*8+1, fItr[1], "Checking 4,0 corner periodic sync")
--      local fItr = field:get(indexer({0,9}))
--      assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
--      local fItr = field:get(indexer({4,9}))
--      assert_equal(4+2*1+1, fItr[1], "Checking 4,9 corner periodic sync")
--   else
--      local fItr = field:get(indexer({3,0}))
--      assert_equal(3+2*8+1, fItr[1], "Checking 3,0 corner periodic sync")
--      local fItr = field:get(indexer({7,0}))
--      assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
--      local fItr = field:get(indexer({3,9}))
--      assert_equal(3+2*1+1, fItr[1], "Checking 3,9 corner periodic sync")
--      local fItr = field:get(indexer({7,9}))
--      assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
--   end

--   -- Check corner ghost cells (cuts={2,1}, periodicDirs={1}).
--   if rank==0 then
--      local fItr = field:get(indexer({0,0}))
--      assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
--      local fItr = field:get(indexer({4,0}))
--      assert_equal(0, fItr[1], "Checking 4,0 corner periodic sync")
--      local fItr = field:get(indexer({0,9}))
--      assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
--      local fItr = field:get(indexer({4,9}))
--      assert_equal(0, fItr[1], "Checking 4,9 corner periodic sync")
--   else
--      local fItr = field:get(indexer({3,0}))
--      assert_equal(0, fItr[1], "Checking 3,0 corner periodic sync")
--      local fItr = field:get(indexer({7,0}))
--      assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
--      local fItr = field:get(indexer({3,9}))
--      assert_equal(0, fItr[1], "Checking 3,9 corner periodic sync")
--      local fItr = field:get(indexer({7,9}))
--      assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
--   end

--   -- Check corner ghost cells (cuts={1,2}, periodicDirs={1}).
--   if rank==0 then
--      local fItr = field:get(indexer({0,0}))
--      assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
--      local fItr = field:get(indexer({7,0}))
--      assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
--      local fItr = field:get(indexer({0,5}))
--      assert_equal(6+2*5+1, fItr[1], "Checking 0,5 corner periodic sync")
--      local fItr = field:get(indexer({7,5}))
--      assert_equal(1+2*5+1, fItr[1], "Checking 7,5 corner periodic sync")
--   else
--      local fItr = field:get(indexer({0,4}))
--      assert_equal(6+2*4+1, fItr[1], "Checking 0,4 corner periodic sync")
--      local fItr = field:get(indexer({7,4}))
--      assert_equal(1+2*4+1, fItr[1], "Checking 7,4 corner periodic sync")
--      local fItr = field:get(indexer({0,9}))
--      assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
--      local fItr = field:get(indexer({7,9}))
--      assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
--   end
end

function test_7(comm)
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 4 then
      log("Not running test_6 as numProcs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
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
end

function test_8(comm)
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   log(string.format("Running shared CartField tests on %d procs", nz))

   local decomp = DecompRegionCalc.CartProd {
      cuts = {1, 1}, useShared = true
   }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }
   local field1 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }   

   field:clear(10.0)
   field1:clear(20.0)

   -- check if loop covers full grid
   local localCount = 0
   for idx in field:localRangeIter() do
      localCount = localCount+1
   end
   local totalCount = allReduceOneInt(localCount)
   assert_equal(grid:localRange():volume(), totalCount, "Checking if total count is correct")
end

function test_9(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 2 then
      log("Not running test_8 as numProcs not exactly 2")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1} }
   
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }
   field:clear(10.25)

   -- write field
   field:write("CartFieldTest_field_p2.bp", 2.5, 50)

   local fieldIn = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
   }
   fieldIn:clear(0.0)

   local tm, fr = fieldIn:read("CartFieldTest_field_p2.bp")

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

-- This test is a repeat of test 3, but communicating device data using CUDA-aware MPI.
function test_10(comm)
   if not GKYL_HAVE_CUDA then return end

   local nz = Mpi.Comm_size(comm)
   if nz ~= 3 then
      log("Not running test_9 as numProcs not exactly 3")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {3} }
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {12},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
      createDeviceCopy = true,
   }

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+1
      fitr[2] = idx[1]+2
      fitr[3] = idx[1]+3
   end

   -- copy stuff to device
   local err = field:copyHostToDevice()
   assert_equal(0, err, "Checking if copy to device worked")

   -- synchronize ghost cells on device
   field:deviceSync() -- sync field

   -- copy data back for checking.
   field:copyDeviceToHost()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+1, fitr[1], "Checking field value")
      assert_equal(idx[1]+2, fitr[2], "Checking field value")
      assert_equal(idx[1]+3, fitr[3], "Checking field value")
   end

   Mpi.Barrier(comm)
end

-- This test is a repeat of test 4, but communicating device data using CUDA-aware MPI.
function test_11(comm)
   if not GKYL_HAVE_CUDA then return end

   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 4 then
      log("Not running test_10 as numProcs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 2},
      syncCorners = true,
      createDeviceCopy = true,
   }

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
      fitr[2] = idx[1]+2*idx[2]+2
      fitr[3] = idx[1]+2*idx[2]+3
   end

   -- copy stuff to device
   local err = field:copyHostToDevice()
   assert_equal(0, err, "Checking if copy to device worked")

   -- synchronize ghost cells on device
   field:deviceSync() -- sync field

   -- copy data back for checking.
   field:copyDeviceToHost()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
      assert_equal(idx[1]+2*idx[2]+2, fitr[2], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
      assert_equal(idx[1]+2*idx[2]+3, fitr[3], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)
end

-- This test is a repeat of test 5, but communicating device data using CUDA-aware MPI.
function test_12(comm)
   if not GKYL_HAVE_CUDA then return end

   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 2 then
      log("Not running test_11 as numProcs not exactly 2")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {500, 500},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1},
      syncCorners = true,
      createDeviceCopy = true,
   }

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
   end

   -- copy stuff to device
   local err = field:copyHostToDevice()
   assert_equal(0, err, "Checking if copy to device worked")

   -- synchronize ghost cells on device
   field:deviceSync() -- sync field

   -- copy data back for checking.
   field:copyDeviceToHost()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)
end

function test_13(comm)
   if not GKYL_HAVE_CUDA then return end

   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 16 then
      log("Not running test_12 as numProcs not exactly 16")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2, 2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0, 0.0},
      upper = {1.0, 1.0, 1.0, 1.0},
      cells = {64, 64, 64, 64},
      decomposition = decomp,
   }
   local field = DataStruct.Field {
      onGrid           = grid,
      numComponents    = 3,
      ghost            = {1, 1},
      syncCorners      = true,
      createDeviceCopy = true,
   }

   local localRange = field:localRange()
   local indexer = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = idx[1]+2*idx[2]+1
   end

   -- Copy stuff to device.
   local err = field:copyHostToDevice()
   assert_equal(0, err, "Checking if copy to device worked")

   -- Synchronize ghost cells on device.
   field:deviceSync() -- Sync field.

   -- Copy data back for checking.
   field:copyDeviceToHost()

   local localExtRange = field:localExtRange():intersect(field:globalRange())
   for idx in localExtRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
   end
   Mpi.Barrier(comm)
end

comm = Mpi.COMM_WORLD
--test_1(comm)
--test_2(comm)
--test_3(comm)
--test_4(comm)
--test_5(comm)
test_6(comm)
--test_7(comm)
--test_8(comm)
--test_9(comm)
--test_10(comm)
--test_11(comm)
--test_12(comm)
--test_13(comm)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
