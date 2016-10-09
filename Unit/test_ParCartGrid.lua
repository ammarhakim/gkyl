-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grid objects
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

function test_1()
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 2 then
      log("Not running test_1 as numProcs not exactly 2")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20},
      decomposition = decomp,
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")

   local localRange = grid:localRange()
   assert_equal(100, localRange:volume(), "Checking sub-domain volume")
end

function test_2()
   local nz = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nz ~= 4 then
      log("Not running test_2 as numProcs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20},
      decomposition = decomp,
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")

   local localRange = grid:localRange()
   assert_equal(50, localRange:volume(), "Checking sub-domain volume")
end

-- Run tests
test_1()
test_2()

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
