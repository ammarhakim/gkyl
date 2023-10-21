-- Gkyl ------------------------------------------------------------------------
--
-- Test for CartField ADIOS IO.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit             = require "Unit"
local Mpi              = require "Comm.Mpi"
local Adios            = require "Io.Adios"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local DataStruct       = require "DataStruct"
local Grid             = require "Grid"
local Lin              = require "Lib.Linalg"
local DecompRegionCalc = require "Lib.CartDecomp"

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then print(msg) end
end

function test_1w(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if nproc > 1 then
      log("test_1w: Not running as number of procs > 1")
      return
   end

   local ad = Adios.init_mpi(comm)

   -- Different grids need a different Adios IO object.
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},  cells = {10, 10},
      upper = {1.0, 1.0},  ioSystem = ad,
   }
   local field = DataStruct.Field {
      onGrid = grid,  numComponents = 1,
      ghost  = {1, 1},
   }
   field:clear(10.5)

   -- I/O object.
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType(),
      metaData = {polyOrder = 0, basisType = "ms",
                  weight = 1.5,}
   }

   adiosIo:write(field, "field.bp", 3.1, 42)

   -- Test writing multiple fields to one file.
   local auxField = DataStruct.Field {
      onGrid = grid,  numComponents = 1,
      ghost = {1, 1},
   }
   auxField:clear(7.5)

   adiosIo:write({field=field,auxField=auxField}, "twoFields.bp", 3.2, 43)

   Adios.finalize(ad)
end

function test_1r(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if nproc > 1 then
      log("test_1r: Not running as number of procs > 1")
      return
   end

   local ad = Adios.init_mpi(comm)

   -- Different grids need a different Adios IO object.
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},  cells = {10, 10},
      upper = {1.0, 1.0},  ioSystem = ad,
   }
   local field = DataStruct.Field {
      onGrid = grid,  numComponents = 1,
      ghost  = {1, 1},
   }

   -- I/O object.
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType()
   }

   local tmStamp, frNum = adiosIo:read(field, "field.bp")

   assert_equal(3.1, tmStamp, "test_1r: Checking field time-stamp")
   assert_equal(42, frNum, "test_1r: Checking field frame number")

   local indexer, localRange = field:indexer(), field:localRange()
   for i = localRange:lower(1),  localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr = field:get(indexer(i, j))
	 assert_equal(10.5, fitr[1], "test_1r: Checking field value")
      end
   end

   -- Test reading multiple fields from one file.
   local auxField = DataStruct.Field {
      onGrid = grid,  numComponents = 1,
      ghost  = {1, 1},
   }

   local tmStamp, frNum = adiosIo:read({field=field,auxField=auxField}, "twoFields.bp")

   assert_equal(3.2, tmStamp, "test_1r: Checking twoFields time-stamp")
   assert_equal(43, frNum, "test_1r: Checking twoFields frame number")

   local indexer, localRange = field:indexer(), field:localRange()
   for i = localRange:lower(1),  localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr    = field:get(indexer(i, j))
	 local auxFitr = auxField:get(indexer(i, j))
	 assert_equal(10.5, fitr[1], "test_1r: Checking twoFields value 1")
	 assert_equal(7.5, auxFitr[1], "test_1r: Checking twoFields value 2")
      end
   end

   Adios.finalize(ad)
end

function test_2w(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   local lower = { 0.0 }
   local upper = { 1.0 }
   local cells = { 20 }

   if math.abs(cells[1]/nproc-math.floor(cells[1]/nproc)) ~= 0 then
      log("test_2w: Not running as numProcs does not evenly divide cells[1]")
      return
   end

   local ad = Adios.init_mpi(comm)

   local decomp = DecompRegionCalc.CartProd { cuts = { nproc } }
   local grid = Grid.RectCart {
      lower = lower,  cells = cells,
      upper = upper,  decomposition = decomp,
      ioSystem = ad,
   }
   local field = DataStruct.Field {
      onGrid = grid,  numComponents = 3,  ghost = { 1, 1 },
   }

   local localRange, indexer = field:localRange(), field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1] = i + 1
      fitr[2] = i + 2
      fitr[3] = i + 3
   end

   -- I/O object.
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType(),
      metaData = {polyOrder = 2, basisType = "ms",
                  weight = 1.5,}
   }
   adiosIo:write(field, "field_decomp.bp", 3.3, 44)

   Adios.finalize(ad)
end

function test_2r(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   local lower = { 0.0 }
   local upper = { 1.0 }
   local cells = { 20 }

   if math.abs(cells[1]/nproc-math.floor(cells[1]/nproc)) ~= 0 then
      log("test_2w: Not running as numProcs does not evenly divide cells[1]")
      return
   end

   local ad = Adios.init_mpi(comm)

   local decomp = DecompRegionCalc.CartProd { cuts = { nproc } }
   local grid = Grid.RectCart {
      lower = lower,  cells = cells,
      upper = upper,  decomposition = decomp,  ioSystem = ad,
   }
   local field = DataStruct.Field {
      onGrid = grid,  numComponents = 3,  ghost = { 1, 1 },
   }
   local fieldIn = DataStruct.Field {
      onGrid = grid,  numComponents = 3,  ghost = { 1, 1 },
   }

   local localRange, indexer = field:localRange(), field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1] = i + 1
      fitr[2] = i + 2
      fitr[3] = i + 3
   end

   -- I/O object.
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType(),
      metaData = {polyOrder = 2, basisType = "ms",
                  weight = 1.5,}
   }
   local tmStamp, frNum = adiosIo:read(fieldIn, "field_decomp.bp")

   assert_equal(3.3, tmStamp, "test_2r: Checking field time-stamp")
   assert_equal(44, frNum, "test_2r: Checking field frame number")

   for i = localRange:lower(1),  localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr = field:get(indexer(i, j))
	 local gitr = fieldIn:get(indexer(i, j))
         for k = 1, field:numComponents() do
            assert_equal(fitr[k], gitr[k], "test_2r: Checking field value")
         end
      end
   end

   Adios.finalize(ad)
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = Lin.IntVec(1), Lin.IntVec(1)
   sendbuf[1] = localv
   Mpi.Allreduce(sendbuf:data(), recvbuf:data(), 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[1]
end

-- Run tests
test_1w(Mpi.COMM_WORLD)
test_1r(Mpi.COMM_WORLD)
test_2w(Mpi.COMM_WORLD)
test_2r(Mpi.COMM_WORLD)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
