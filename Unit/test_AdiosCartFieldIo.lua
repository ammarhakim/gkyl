-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS IO
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
   local ad = Adios.init_mpi(comm)

   -- Different grids need a different Adios IO object.
   local ad_io = Adios.declare_io(ad, "gridio")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},  cells = {10, 10},
      upper = {1.0, 1.0},  ioSystem = ad_io,
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
   local ad = Adios.init_mpi(comm)

   -- Different grids need a different Adios IO object.
   local ad_io = Adios.declare_io(ad, "gridio")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},  cells = {10, 10},
      upper = {1.0, 1.0},  ioSystem = ad_io,
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

function allReduceOneInt(localv)
   local sendbuf, recvbuf = Lin.IntVec(1), Lin.IntVec(1)
   sendbuf[1] = localv
   Mpi.Allreduce(sendbuf:data(), recvbuf:data(), 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[1]
end

-- Run tests
test_1w(Mpi.COMM_WORLD)
test_1r(Mpi.COMM_WORLD)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
