-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition of Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local Mpi = require "Comm.Mpi"

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
      log("Test for Split_comm not run as number of procs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   assert_equal(false, decomp:isShared(), "Checking if decomp is shared")

   local worldComm = decomp:comm()
   local shmComm = decomp:sharedComm()
   local nodeComm = decomp:nodeComm()

   assert_equal(4, Mpi.Comm_size(worldComm), "Checking size for global comm")
   assert_equal(1, Mpi.Comm_size(shmComm), "Checking size for share comm")
   if Mpi.Is_comm_valid(nodeComm) then
      assert_equal(4, Mpi.Comm_size(nodeComm), "Checking size of node comm")
   end

   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))
   assert_equal(4, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- fetch domains and do sanity checks
   local v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):volume()
   end
   assert_equal(100, v, "Checking volume of decomposedRgn")
end

function test_2(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for Split_comm not run as number of procs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1}, useShared = true }
   assert_equal(true, decomp:isShared(), "Checking if decomp is shared")

   local worldComm = decomp:comm()
   local shmComm = decomp:sharedComm()
   local nodeComm = decomp:nodeComm()

   assert_equal(4, Mpi.Comm_size(worldComm), "Checking size for global comm")
   assert_equal(4, Mpi.Comm_size(shmComm), "Checking size for share comm")
   if Mpi.Is_comm_valid(nodeComm) then
      assert_equal(1, Mpi.Comm_size(nodeComm), "Checking size of node comm")
   end

   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))
   assert_equal(1, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- fetch domains and do sanity checks
   local v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):volume()
   end
   assert_equal(100, v, "Checking volume of decomposedRgn")   
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
