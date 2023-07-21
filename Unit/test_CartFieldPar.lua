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

function copyFieldGlobalFromLocal(fieldGlobal, fieldGlobalBuffer, fieldLocal, fieldLocalBuffer)
   local function copyField(f1, f2) -- copies contents of f2 into f1
      -- To copy field Global to field Local, input f1-FieldGlobal, f2-FieldLocal
      f1:copyRangeToRange(f2, f1:localRange(), f2:localRange())
   end
   copyField(fieldLocalBuffer, fieldLocal)
   -- Perform the Mpi.Allgather operation
   Mpi.Allgather(fieldLocalBuffer:dataPointer(), fieldLocalBuffer:size(), fieldLocalBuffer:elemCommType(),
      fieldGlobalBuffer:dataPointer(), fieldLocalBuffer:size(), fieldGlobalBuffer:elemCommType(), comm)
   copyField(fieldGlobal, fieldGlobalBuffer)
   return
end

function copyFieldLocalFromGlobal(fieldLocal, fieldGlobal)
   fieldLocal:copyRangeToRange(fieldGlobal, fieldLocal:localRange(), fieldLocal:localRange())
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

   if GKYL_USE_GPU then
      field:copyHostToDevice()
   end

   field:sync()

--   field._zero:clear(0)

   if GKYL_USE_GPU then
      field:copyDeviceToHost()
   end

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

   if GKYL_USE_GPU then
      field:copyHostToDevice()
   end

   field:sync()

--   field._zero:clear(0)

   if GKYL_USE_GPU then
      field:copyDeviceToHost()
   end

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
      syncCorners   = false,
   }

   local localRange = field:localRange()
   local indexer    = field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      fitr[1]    = idx[1]+2*idx[2]+1
   end

   if GKYL_USE_GPU then
      field:copyHostToDevice()
   end

   field:sync()

--   field._zero:clear(0)

   if GKYL_USE_GPU then
      field:copyDeviceToHost()
   end

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
      log("Not running test_6 as numProcs not exactly 2")
      return
   end

   local testCuts = {{2,1},{1,2},{2,1},{2,1},{1,2}}
   local testPeriodicDirs = {{1,2},{1,2},{2},{1},{1}} 

   local decomp, grid, field = {}, {}, {}
   for tI = 1, #testCuts do
      decomp[tI] = DecompRegionCalc.CartProd { cuts = testCuts[tI] }

      grid[tI] = Grid.RectCart {
         lower = {0.0, 0.0},
         upper = {1.0, 1.0},
         cells = {6, 8},
         decomposition = decomp[tI],
         periodicDirs  = testPeriodicDirs[tI],
      }
      field[tI] = DataStruct.Field {
         onGrid        = grid[tI],
         numComponents = 1,
         ghost         = {1, 1},
         syncCorners   = true,
      }

      local localRange = field[tI]:localRange()
      local indexer    = field[tI]:genIndexer()
      for idx in localRange:colMajorIter() do
         local fitr = field[tI]:get(indexer(idx))
         fitr[1]    = idx[1]+2*idx[2]+1
      end

      if GKYL_USE_GPU then
         field[tI]:copyHostToDevice()
      end

      field[tI]:sync()

      if GKYL_USE_GPU then
         field[tI]:copyDeviceToHost()
      end

      local localExtRange = field[tI]:localExtRange():intersect(field[tI]:globalRange())
      for idx in localExtRange:colMajorIter() do
         local fitr = field[tI]:get(indexer(idx))
         assert_equal(idx[1]+2*idx[2]+1, fitr[1], string.format("Checking field value at (%d, %d)", idx[1], idx[2]))
      end
      Mpi.Barrier(comm)

      local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

      if tI == 1 then
         -- Check corner ghost cells (cuts={2,1}, periodicDirs={1,2}).
         if rank==0 then
            local fItr = field[tI]:get(indexer({0,0}))
            assert_equal(6+2*8+1, fItr[1], "Checking 0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,0}))
            assert_equal(4+2*8+1, fItr[1], "Checking 4,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9}))
            assert_equal(6+2*1+1, fItr[1], "Checking 0,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,9}))
            assert_equal(4+2*1+1, fItr[1], "Checking 4,9 corner periodic sync")
         else
            local fItr = field[tI]:get(indexer({3,0}))
            assert_equal(3+2*8+1, fItr[1], "Checking 3,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0}))
            assert_equal(1+2*8+1, fItr[1], "Checking 7,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({3,9}))
            assert_equal(3+2*1+1, fItr[1], "Checking 3,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9}))
            assert_equal(1+2*1+1, fItr[1], "Checking 7,9 corner periodic sync")
         end
      
      elseif tI == 2 then
         -- Check corner ghost cells (cuts={1,2}, periodicDirs={1,2}).
         if rank==0 then
            local fItr = field[tI]:get(indexer({0,0}))
            assert_equal(6+2*8+1, fItr[1], "Checking 0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0}))
            assert_equal(1+2*8+1, fItr[1], "Checking 7,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,5}))
            assert_equal(6+2*5+1, fItr[1], "Checking 0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,5}))
            assert_equal(1+2*5+1, fItr[1], "Checking 7,5 corner periodic sync")
         else
            local fItr = field[tI]:get(indexer({0,4}))
            assert_equal(6+2*4+1, fItr[1], "Checking 0,4 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,4}))
            assert_equal(1+2*4+1, fItr[1], "Checking 7,4 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9}))
            assert_equal(6+2*1+1, fItr[1], "Checking 0,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9}))
            assert_equal(1+2*1+1, fItr[1], "Checking 7,9 corner periodic sync")
         end
      
      elseif tI == 3 then
         -- Check corner ghost cells (cuts={2,1}, periodicDirs={2}).
         if rank==0 then
            local fItr = field[tI]:get(indexer({0,0}))
            assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,0}))
            assert_equal(4+2*8+1, fItr[1], "Checking 4,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9}))
            assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,9}))
            assert_equal(4+2*1+1, fItr[1], "Checking 4,9 corner periodic sync")
         else
            local fItr = field[tI]:get(indexer({3,0}))
            assert_equal(3+2*8+1, fItr[1], "Checking 3,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0}))
            assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({3,9}))
            assert_equal(3+2*1+1, fItr[1], "Checking 3,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9}))
            assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
         end
      
      elseif tI == 4 then
         -- Check corner ghost cells (cuts={2,1}, periodicDirs={1}).
         if rank==0 then
            local fItr = field[tI]:get(indexer({0,0}))
            assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,0}))
            assert_equal(0, fItr[1], "Checking 4,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9}))
            assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,9}))
            assert_equal(0, fItr[1], "Checking 4,9 corner periodic sync")
         else
            local fItr = field[tI]:get(indexer({3,0}))
            assert_equal(0, fItr[1], "Checking 3,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0}))
            assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({3,9}))
            assert_equal(0, fItr[1], "Checking 3,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9}))
            assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
         end
      
      elseif tI == 5 then
         -- Check corner ghost cells (cuts={1,2}, periodicDirs={1}).
         if rank==0 then
            local fItr = field[tI]:get(indexer({0,0}))
            assert_equal(0, fItr[1], "Checking 0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0}))
            assert_equal(0, fItr[1], "Checking 7,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,5}))
            assert_equal(6+2*5+1, fItr[1], "Checking 0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,5}))
            assert_equal(1+2*5+1, fItr[1], "Checking 7,5 corner periodic sync")
         else
            local fItr = field[tI]:get(indexer({0,4}))
            assert_equal(6+2*4+1, fItr[1], "Checking 0,4 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,4}))
            assert_equal(1+2*4+1, fItr[1], "Checking 7,4 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9}))
            assert_equal(0, fItr[1], "Checking 0,9 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9}))
            assert_equal(0, fItr[1], "Checking 7,9 corner periodic sync")
         end
      end
   end
end

function test_7(comm)
   -- Check sync in 3D, and the corners shared by the two MPI ranks.
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 2 then
      log("Not running test_7 as numProcs not exactly 2")
      return
   end

   local testCuts = {{2,1,1}}
   local testPeriodicDirs = {{1,2,3}} 

   local decomp, grid, field = {}, {}, {}
   for tI = 1, #testCuts do
      decomp[tI] = DecompRegionCalc.CartProd { cuts = testCuts[tI] }

      grid[tI] = Grid.RectCart {
         lower = {0.0, 0.0, 0.0},
         upper = {1.0, 1.0, 1.0},
         cells = {6, 8, 4},
         decomposition = decomp[tI],
         periodicDirs  = testPeriodicDirs[tI],
      }
      field[tI] = DataStruct.Field {
         onGrid        = grid[tI],
         numComponents = 1,
         ghost         = {1, 1},
         syncCorners   = true,
      }

      local localRange = field[tI]:localRange()
      local indexer    = field[tI]:genIndexer()
      for idx in localRange:colMajorIter() do
         local fitr = field[tI]:get(indexer(idx))
         fitr[1]    = idx[1]+2*idx[2]+3*idx[3]+1
      end

      if GKYL_USE_GPU then
         field[tI]:copyHostToDevice()
      end

      field[tI]:sync()

      if GKYL_USE_GPU then
         field[tI]:copyDeviceToHost()
      end

      local localExtRange = field[tI]:localExtRange():intersect(field[tI]:globalRange())
      for idx in localExtRange:colMajorIter() do
         local fitr = field[tI]:get(indexer(idx))
         assert_equal(idx[1]+2*idx[2]+3*idx[3]+1, fitr[1], string.format("Checking field value at (%d, %d, %d)", idx[1], idx[2], idx[3]))
      end
      Mpi.Barrier(comm)

      local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

      if tI == 1 then
         -- Check "corner" ghost cells (cuts={2,1}, periodicDirs={1,2}).
         if rank==0 then
            -- Check corners.
            local fItr = field[tI]:get(indexer({0,0,0}))
            assert_equal(6+2*8+3*4+1, fItr[1], "Checking 0,0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,0,0}))
            assert_equal(4+2*8+3*4+1, fItr[1], "Checking 4,0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9,0}))
            assert_equal(6+2*1+3*4+1, fItr[1], "Checking 0,9,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,9,0}))
            assert_equal(4+2*1+3*4+1, fItr[1], "Checking 4,9,0 corner periodic sync")

            local fItr = field[tI]:get(indexer({0,0,5}))
            assert_equal(6+2*8+3*1+1, fItr[1], "Checking 0,0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,0,5}))
            assert_equal(4+2*8+3*1+1, fItr[1], "Checking 4,0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({0,9,5}))
            assert_equal(6+2*1+3*1+1, fItr[1], "Checking 0,9,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({4,9,5}))
            assert_equal(4+2*1+3*1+1, fItr[1], "Checking 4,9,5 corner periodic sync")

            -- Check ghost edges.
            for i=1,3 do
               local fItr = field[tI]:get(indexer({i,0,0}))
               assert_equal(i+2*8+3*4+1, fItr[1], string.format("Checking %d,0,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,9,0}))
               assert_equal(i+2*1+3*4+1, fItr[1], string.format("Checking %d,9,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,0,5}))
               assert_equal(i+2*8+3*1+1, fItr[1], string.format("Checking %d,0,5 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,9,5}))
               assert_equal(i+2*1+3*1+1, fItr[1], string.format("Checking %d,9,5 corner periodic sync",i))
            end
            for i=1,8 do
               local fItr = field[tI]:get(indexer({0,i,0}))
               assert_equal(6+2*i+3*4+1, fItr[1], string.format("Checking 0,%d,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({0,i,5}))
               assert_equal(6+2*i+3*1+1, fItr[1], string.format("Checking 0,%d,5 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({4,i,0}))
               assert_equal(4+2*i+3*4+1, fItr[1], string.format("Checking 4,%d,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({4,i,5}))
               assert_equal(4+2*i+3*1+1, fItr[1], string.format("Checking 4,%d,5 corner periodic sync",i))
            end
            for i=1,4 do
               local fItr = field[tI]:get(indexer({0,0,i}))
               assert_equal(6+2*8+3*i+1, fItr[1], string.format("Checking 0,0,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({4,0,i}))
               assert_equal(4+2*8+3*i+1, fItr[1], string.format("Checking 4,0,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({0,9,i}))
               assert_equal(6+2*1+3*i+1, fItr[1], string.format("Checking 0,9,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({4,9,i}))
               assert_equal(4+2*1+3*i+1, fItr[1], string.format("Checking 4,9,%d corner periodic sync",i))
            end
         else
            -- Check corners.
            local fItr = field[tI]:get(indexer({3,0,0}))
            assert_equal(3+2*8+3*4+1, fItr[1], "Checking 3,0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0,0}))
            assert_equal(1+2*8+3*4+1, fItr[1], "Checking 7,0,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({3,9,0}))
            assert_equal(3+2*1+3*4+1, fItr[1], "Checking 3,9,0 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9,0}))
            assert_equal(1+2*1+3*4+1, fItr[1], "Checking 7,9,0 corner periodic sync")

            local fItr = field[tI]:get(indexer({3,0,5}))
            assert_equal(3+2*8+3*1+1, fItr[1], "Checking 3,0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,0,5}))
            assert_equal(1+2*8+3*1+1, fItr[1], "Checking 7,0,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({3,9,5}))
            assert_equal(3+2*1+3*1+1, fItr[1], "Checking 3,9,5 corner periodic sync")
            local fItr = field[tI]:get(indexer({7,9,5}))
            assert_equal(1+2*1+3*1+1, fItr[1], "Checking 7,9,5 corner periodic sync")

            -- Check ghost edges.
            for i=4,6 do
               local fItr = field[tI]:get(indexer({i,0,0}))
               assert_equal(i+2*8+3*4+1, fItr[1], string.format("Checking %d,0,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,9,0}))
               assert_equal(i+2*1+3*4+1, fItr[1], string.format("Checking %d,9,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,0,5}))
               assert_equal(i+2*8+3*1+1, fItr[1], string.format("Checking %d,0,5 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({i,9,5}))
               assert_equal(i+2*1+3*1+1, fItr[1], string.format("Checking %d,9,5 corner periodic sync",i))
            end
            for i=1,8 do
               local fItr = field[tI]:get(indexer({3,i,0}))
               assert_equal(3+2*i+3*4+1, fItr[1], string.format("Checking 3,%d,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({3,i,5}))
               assert_equal(3+2*i+3*1+1, fItr[1], string.format("Checking 3,%d,5 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({7,i,0}))
               assert_equal(1+2*i+3*4+1, fItr[1], string.format("Checking 7,%d,0 corner periodic sync",i))
               local fItr = field[tI]:get(indexer({7,i,5}))
               assert_equal(1+2*i+3*1+1, fItr[1], string.format("Checking 7,%d,5 corner periodic sync",i))
            end
            for i=1,4 do
               local fItr = field[tI]:get(indexer({3,0,i}))
               assert_equal(3+2*8+3*i+1, fItr[1], string.format("Checking 3,0,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({7,0,i}))
               assert_equal(1+2*8+3*i+1, fItr[1], string.format("Checking 7,0,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({3,9,i}))
               assert_equal(3+2*1+3*i+1, fItr[1], string.format("Checking 3,9,%d corner periodic sync",i))
               local fItr = field[tI]:get(indexer({7,9,i}))
               assert_equal(1+2*1+3*i+1, fItr[1], string.format("Checking 7,9,%d corner periodic sync",i))
            end
         end
      end
   end
end
      
function test_8(comm)
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 4 then
      log("Not running test_8 as numProcs not exactly 4")
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

   if GKYL_USE_GPU then
      field1:copyHostToDevice()
   end

   -- accumulate stuff
   field:accumulate(1.0, field1, 2.0, field1)

   if GKYL_USE_GPU then
      field:copyDeviceToHost()
   end

   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      assert_equal(10+3*(idx[1]+2*idx[2]+1), fitr[1], "Checking field value")
      assert_equal(10+3*(idx[1]+2*idx[2]+2), fitr[2], "Checking field value")
      assert_equal(10+3*(idx[1]+2*idx[2]+3), fitr[3], "Checking field value")
   end   
end

function test_9(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 2 then
      log("Not running test_9 as numProcs not exactly 2")
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

function test_10(comm)
   -- Test sending/receving a whole CartField with Isend/Irecv and Wait.
   local defaultCommSize = 4
   local nz = Mpi.Comm_size(comm)
   if nz ~= defaultCommSize then
      log(string.format("Not running test_10 as numProcs not exactly %d",defaultCommSize))
      return
   end

   local worldRank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   local numCutsConf = 2 --defaultCommSize
   
   local confColor = math.floor(worldRank/numCutsConf)
   local confComm  = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, worldRank)
   local speciesColor = worldRank % numCutsConf
   local speciesComm  = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, worldRank)
   local confRank = Mpi.Comm_rank(confComm)
   local speciesRank = Mpi.Comm_rank(speciesComm)

   local rank = numCutsConf==2 and speciesRank or confRank
   local comm = numCutsConf==2 and speciesComm or confComm

   local decomp = DecompRegionCalc.CartProd { cuts = {numCutsConf}, comm = confComm }
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {12},
      decomposition = decomp,
   }
   local field0 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   local field1 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   if rank==0 then field0:clear(0.3) end
   if rank==1 then field1:clear(1.5) end

   local recvReq, sendReq = Mpi.Request(), Mpi.Request()
   local recvStat, sendStat = Mpi.Status(), Mpi.Status()
   local tag = 32
   if rank==0 then
      local _ = Mpi.Irecv(field1:dataPointer(), field1:size(), field1:elemCommType(),
                          1, tag, comm, recvReq)
      local _ = Mpi.Isend(field0:dataPointer(), field0:size(), field0:elemCommType(),
                          1, tag, comm, sendReq)
   else
      local _ = Mpi.Irecv(field0:dataPointer(), field0:size(), field0:elemCommType(),
                          0, tag, comm, recvReq)
      local _ = Mpi.Isend(field1:dataPointer(), field1:size(), field1:elemCommType(),
                          0, tag, comm, sendReq)
   end

   local _ = Mpi.Wait(recvReq, recvStat)

   local field = rank==0 and field1 or field0
   local value = rank==0 and 1.5 or 0.3
   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      assert_equal(value, fitr[1], "test_10: Checking field value")
      assert_equal(value, fitr[2], "test_10: Checking field value")
      assert_equal(value, fitr[3], "test_10: Checking field value")
   end

   local _ = Mpi.Wait(sendReq, sendStat)
   Mpi.Barrier(comm)
end

function test_11(comm)
   --Test the establishment of a global grid and field
   local nz = Mpi.Comm_size(comm)
   if nz ~= 2 then
      log("Not running test_11 as numProcs not exactly 2")
      return
   end
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

   local decomp = DecompRegionCalc.CartProd { cuts = {2} }   
   local noDecomp = DecompRegionCalc.CartProd { 
      cuts = {1},
      __serTesting = true, --try without this flag
    }
   local grid = Grid.RectCart {--located in Grid/RectCart.lua
      lower = {0.0},
      upper = {1.0},
      cells = {10},
      decomposition = decomp,-- what does decomp do?
   }
   local field = DataStruct.Field { --located in DataStruct/CartField.lua
      onGrid        = grid,
      numComponents = 3,--what does 3 mean? 3 what? 3 DG polynomials?
      ghost         = {1, 1},
   }
   local gridGlobal = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
      decomposition = noDecomp,
   }
   local fieldGlobal = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = 3,
      ghost         = {1, 1},
   }

   --Local field checks
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

   --Global field checks
   assert_equal(1, fieldGlobal:ndim(), "Checking dimensions")
   assert_equal(1, fieldGlobal:lowerGhost(), "Checking global lower ghost")
   assert_equal(1, fieldGlobal:upperGhost(), "Checking global upper ghost") -- What should this value be?
   assert_equal(fieldGlobal:localExtRange():volume()*3, fieldGlobal:size(), "Checking global size")

   assert_equal(fieldGlobal:defaultLayout(), fieldGlobal:layout(), "Checking global layout")

   local localRangeGlobal = fieldGlobal:localRange()
   assert_equal(1, localRangeGlobal:lower(1), "Checking global range lower")
   assert_equal(10, localRangeGlobal:upper(1), "Checking global range upper")

   local localExtRangeGlobal = fieldGlobal:localExtRange()
   assert_equal(0, localExtRangeGlobal:lower(1), "Checking global ext range lower")
   assert_equal(11, localExtRangeGlobal:upper(1), "Checking global ext range upper")

   local indexer = fieldGlobal:indexer() -- there is an error here
   for i = localRangeGlobal:lower(1), localRangeGlobal:upper(1) do
      local fitr = fieldGlobal:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end

   for i = localRangeGlobal:lower(1), localRangeGlobal:upper(1) do
      local fitr = fieldGlobal:get(indexer(i))
      assert_equal(i+1, fitr[1], "Checking global field value")
      assert_equal(i+2, fitr[2], "Checking global field value")
      assert_equal(i+3, fitr[3], "Checking global field value")
   end
   Mpi.Barrier(comm)
end

function test_12(comm)
   local nz = Mpi.Comm_size(comm)
   if nz ~= 2 then
      log("Not running test_12 as numProcs not exactly 2")
      return
   end
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

   local decomp = DecompRegionCalc.CartProd { cuts = {2} }   
   local noDecomp = DecompRegionCalc.CartProd { 
      cuts = {1},
      __serTesting = true,
    }
   local grid = Grid.RectCart {--located in Grid/RectCart.lua
      lower = {0.0},
      upper = {1.0},
      cells = {10},
      decomposition = decomp,-- what does decomp do?
   }
   local gridGlobal = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
      decomposition = noDecomp,
   }
   assert_equal(1, gridGlobal:globalRange():lower(1), "Checking range lower")
   assert_equal(10, gridGlobal:globalRange():upper(1), "Checking range upper")
   assert_equal(1,  gridGlobal:localRange():lower(1), "Checking range lower")
   assert_equal(10, gridGlobal:localRange():upper(1), "Checking range upper")
   Mpi.Barrier(comm)
end

function test_copyFieldGlobalFromLocal_1D(comm)
   -- Test Allgathering fieldLocals into a fieldGlobal
   -- works in 1D case

   local sz = Mpi.Comm_size(comm)
   if sz < 2 then
      log("Not running test_13 as numProcs less than 2")
      return
   end

   local rank = Mpi.Comm_rank(comm)
   local nGhost = 1
   local nComponent = 3

   local nCells = {sz}
   local lowerBound = {0.0}
   local upperBound = {0.0}

   -- define grids and decomposition cuts
   local decomp = DecompRegionCalc.CartProd { 
	  cuts = {sz}
   }   
   local noDecomp = DecompRegionCalc.CartProd { 
      cuts = {1},
      __serTesting = true,
    }
   local gridLocal = Grid.RectCart {--located in Grid/RectCart.lua
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = decomp,-- what does decomp do?
   }
   local gridGlobal = Grid.RectCart {
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = noDecomp,
   }

   local fieldLocal = DataStruct.Field {
      onGrid        = gridLocal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }
   local fieldGlobal = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }

   -- Buffer fields have no Ghost cells for use in MPI Allreduce
   local fieldLocalBuffer = DataStruct.Field {
      onGrid        = gridLocal,
      numComponents = nComponent,
      ghost         = {0, 0},
   }
   local fieldGlobalBuffer = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = nComponent,
      ghost         = {0, 0},
   }

   -- set local field to unique value for each processor
   fieldLocal:clear(rank)

   copyFieldGlobalFromLocal(fieldGlobal, fieldGlobalBuffer, fieldLocal, fieldLocalBuffer)

   -- Verify the correctness of the gathered data
   local field = fieldGlobal

   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local value = i-1
      local fitr = field:get(indexer(i))
      assert_equal(value, fitr[1], "test_10: Checking field value")
      assert_equal(value, fitr[2], "test_10: Checking field value")
      assert_equal(value, fitr[3], "test_10: Checking field value")

   end

   -- Barrier synchronization to ensure all processes complete the test
   Mpi.Barrier(comm)
end

function test_copyFieldGlobalFromLocal_2D(comm)
   -- Test Allgathering fieldLocals into a fieldGlobal
   -- for multiDimension

   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Not running test_14 as numProcs not exactly 4")
      return
   end

   local rank = Mpi.Comm_rank(comm)
   local nGhost = 1
   local nComponent = 3

   local nCells = {10,10}
   local lowerBound = {0.0, 0.0}
   local upperBound = {1.0, 1.0}

   -- define grids and decomposition cuts
   local decomp = DecompRegionCalc.CartProd {
	  cuts = {2,2}
   }
   local noDecomp = DecompRegionCalc.CartProd {
      cuts = {1,1},
      __serTesting = true, 
    }

   local gridLocal = Grid.RectCart {
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = decomp,-- what does decomp do?
   }
   local gridGlobal = Grid.RectCart {
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = noDecomp,
   }

   -- define a local field and global field on each process
   local fieldLocal = DataStruct.Field {
      onGrid        = gridLocal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }
   local fieldGlobal = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }

   -- Buffer fields have no Ghost cells for use in MPI Allreduce
   local fieldLocalBuffer = DataStruct.Field {
      onGrid        = gridLocal,
      numComponents = nComponent,
      ghost         = {0, 0},
   }
   local fieldGlobalBuffer = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = nComponent,
      ghost         = {0, 0},
   }

   -- set local field to unique value for each processor
   fieldLocal:clear(rank)

   -- copy the localField to the globalField
   copyFieldGlobalFromLocal(fieldGlobal, fieldGlobalBuffer, fieldLocal, fieldLocalBuffer)

   -- Verify the correctness of the gathered data
   local field = fieldGlobal

   local indexer = field:genIndexer()
   for i in field:localRangeIter() do
      local idx = indexer(i) - 14
      local fitr = field:get(indexer(i))
      local j = 2
      if idx < 25 then
         value = 0
         assert_equal(value, fitr[j], "test_14: Checking field value")
      elseif idx < 50 then
         value = 1
         assert_equal(value, fitr[j], "test_14: Checking field value")
      elseif idx < 75 then
         value = 2
         assert_equal(value, fitr[j], "test_14: Checking field value")
      elseif idx < 100 then
         value = 3
         assert_equal(value, fitr[j], "test_14: Checking field value")
      else
         value = 0
         assert_equal(value, fitr[j], "test_14: Checking field value")
      end
   end
   -- Barrier synchronization to ensure all processes complete the test
   Mpi.Barrier(comm)
end

-- Test for the copyField function
function test_copyFieldLocalFromGlobal()
   -- Initialize the global field and grid
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Not running test_14 as numProcs not exactly 4")
      return
   end

   local rank = Mpi.Comm_rank(comm)
   local nGhost = 1
   local nComponent = 1

   local nCells = {2,2}
   local lowerBound = {0.0, 0.0}
   local upperBound = {1.0, 1.0}

   -- define grids and decomposition cuts
   local decomp = DecompRegionCalc.CartProd {
	  cuts = {2,2}
   }
   local noDecomp = DecompRegionCalc.CartProd {
      cuts = {1,1},
      __serTesting = true,
    }

   local gridLocal = Grid.RectCart {
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = decomp,-- what does decomp do?
   }
   local gridGlobal = Grid.RectCart {
      lower = lowerBound,
      upper = upperBound,
      cells = nCells,
      decomposition = noDecomp,
   }

   -- define a local field and global field on each process
   local fieldLocal = DataStruct.Field {
      onGrid        = gridLocal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }
   local fieldGlobal = DataStruct.Field {
      onGrid        = gridGlobal,
      numComponents = nComponent,
      ghost         = {nGhost, nGhost},
   }
   function fPrint(field,msg)
      local indexer = field:genIndexer()
      for i in field:localRangeIter() do
         local fitr = field:get(indexer(i))
         print(msg,": rank, val", rank, fitr[2])
      end
   end


   --Set the fieldGlobal to 5, but the local field to something else
   fieldGlobal:clear(5.0)
   fieldLocal:clear(3.0)
   
   local testRank = 1
   if rank == testRank then fPrint(fieldLocal,"fieldLocal")   end
   if rank == testRank then fPrint(fieldGlobal,"fieldGlobal") end
--if rank == testRank then print(localRangePrint:lower(1), localRangePrint:lower(2)) end
--if rank == testRank then print(localRangePrint:upper(1), localRangePrint:upper(2)) end


   if rank == testRank then print("I am copying") end
   -- Copy the global field to the local field using copyField()
   copyFieldLocalFromGlobal(fieldLocal, fieldGlobal)
   Mpi.Barrier(comm)
   if rank == testRank then fPrint(fieldLocal,"fieldLocal")   end
   if rank == testRank then fPrint(fieldGlobal,"fieldGlobal") end

   Mpi.Barrier(comm)

   -- Verify that the local field now has the same data as the global field
   local indexer = fieldLocal:genIndexer()
   for idx in fieldLocal:localRangeIter() do
      local fitrIn = fieldLocal:get(indexer(idx))
      for k = 1, fieldLocal:numComponents() do
	      assert_equal(5, fitrIn[k], "Checking if field read correctly")
      end
   end
   Mpi.Barrier(comm)
end



comm = Mpi.COMM_WORLD
-- test_1(comm)
-- test_2(comm)
-- test_3(comm)
-- test_4(comm)
-- test_5(comm)
-- test_6(comm)
-- test_7(comm)
-- test_8(comm)
-- test_9(comm)
-- test_10(comm)
-- test_11(comm)
-- test_12(comm)
--test_copyFieldGlobalFromLocal_1D(comm)
--test_copyFieldGlobalFromLocal_2D(comm)
test_copyFieldLocalFromGlobal(comm)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
