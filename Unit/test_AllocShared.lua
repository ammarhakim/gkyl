-- Gkyl ------------------------------------------------------------------------
--
-- Test for shared memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Mpi = require "Comm.Mpi"
local AllocShared = require "Lib.AllocShared"

local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function allReduceOneInt(comm, localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, comm)
   return recvbuf[0]
end

function makeLogger(comm, logRank)
   local rank = Mpi.Comm_rank(comm)
   return function (msg)
      if rank == logRank then 
	 io.write(msg); io.write("\n")
      end
   end
end
log = makeLogger(Mpi.COMM_WORLD, 0)

function test_1(comm)
   local vec = AllocShared.Double(comm, 22)

   assert_equal(true, vec:elemType() == typeof("double"), "Testing elem type")
   assert_equal(sizeof("double"), vec:elemSize(), "Testing elem size")
   assert_equal(22, vec:size(), "Checking total size")

   local glbSize = allReduceOneInt(comm, vec:localSize())
   assert_equal(vec:size(), glbSize, "Checking sum of local sizes")

   for i = vec:lower(), vec:upper() do
      vec[i] = i
   end

   local rnk = Mpi.Comm_rank(comm)
   if rnk == 0 then
      for i = 1, vec:size() do
	 assert_equal(i, vec[i], "Testing if parallel for worked")
      end
   end
end

-- run tests
worldSz = Mpi.Comm_size(Mpi.COMM_WORLD)
shmComm = Mpi.Comm_split_type(Mpi.COMM_WORLD, Mpi.COMM_TYPE_SHARED, 0, Mpi.INFO_NULL)
shmSz = Mpi.Comm_size(shmComm)
log(string.format("Total number of procs %d; per-node procs %d", worldSz, shmSz))

test_1(shmComm)

totalFail = allReduceOneInt(Mpi.COMM_WORLD, stats.fail)
totalPass = allReduceOneInt(Mpi.COMM_WORLD, stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
