-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grids in parallel
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local DecompRegionCalc = require "Lib.CartDecomp"

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

function test_1(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {1.0, 1.0},
      cells = {100, 200},
      decomposition = decomp,
   }

   local localRange = grid:localRange()
   assert_equal(50*100, localRange:volume(), "Checking volume of local ranges")
   assert_equal(rnk+1, grid:subGridId(), "Checking sub-grid ID is correct")
   
end

function test_2(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1}, useShared = true }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {1.0, 1.0},
      cells = {100, 200},
      decomposition = decomp,
   }

   local localRange = grid:localRange()
   assert_equal(100*200, localRange:volume(), "Checking volume of local ranges")

   assert_equal(1, grid:subGridId(), "Checking sub-grid ID is correct")
end

-- Run tests
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
