-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS IO
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Unit = require "Unit"

local ffi  = require "ffi"
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

function test_1w(comm)
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1},
   }
   field:clear(10.5)

   -- I/O object
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType(),
      metaData = {
	 polyOrder = 2, basisType = "ms",
	 weight = 1.5,
      }
   }

   adiosIo:write(field, "field.bp", 3.1, 42)
end

function test_1r(comm)
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1},
   }

   -- I/O object
   local adiosIo = AdiosCartFieldIo {
      elemType = field:elemType()
   }

   local tmStamp, frNum = adiosIo:read(field, "field.bp")

   assert_equal(3.1, tmStamp, "Checking time-stamp")
   assert_equal(42, frNum, "Checking frame number")

   local indexer = field:indexer()
   local localRange = field:localRange()
   for i = localRange:lower(1),  localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 local fitr = field:get(indexer(i, j))
	 assert_equal(10.5, fitr[1], "Checking field value")
      end
   end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
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
