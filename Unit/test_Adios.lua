-- Gkyl ------------------------------------------------------------------------
--
-- Test for MPI wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"

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
   local nproc = Mpi.Comm_size(comm)
   if nproc > 1 then
      log("Not running test_1 as number of procs > 1")
   end
   
   local rank = Mpi.Comm_rank(comm)

   Adios.init_noxml(comm)
   Adios.set_max_buffer_size(16) -- 16 MB chunks

   -- create group and set I/O method
   local grpId = Adios.declare_group("Moments", "", Adios.flag_no)
   Adios.select_method(grpId, "MPI", "", "")

   -- define attributes
   local cells = new("int[2]")
   cells[0] = 10; cells[1] = 20
   Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, 2, cells)

   local lower = new("double[2]")
   lower[0] = 0.0; lower[1] = 0.0
   Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, 2, lower)

   local upper = new("double[2]")
   upper[0] = 1.5; upper[1] = 1.5
   Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, 2, upper)
   
   -- define variables
   Adios.define_var(grpId, "NX", "", Adios.integer, "", "", "")
   Adios.define_var(grpId, "NY", "", Adios.integer, "", "", "")
   Adios.define_var(grpId, "temperature", "", Adios.double, "NX", "", "")
   -- open file to write out group
   local fd = Adios.open("Moments", "adios-test-1.bp", "w", comm)

   -- compute group size
   local nx = new("int[1]"); nx[0] = 100
   local ny = new("int[1]"); ny[0] = 20
   local temperature = new("double[?]", nx[0])

   for i = 0, nx[0]-1 do
      temperature[i] = i
   end
   
   local grpSize = 2*sizeof("int") + nx[0]*sizeof("double")
   local totalSize = Adios.group_size(fd, grpSize)

   -- write data
   Adios.write(fd, "NX", nx)
   Adios.write(fd, "NY", ny)
   Adios.write(fd, "temperature", temperature)
   
   Adios.close(fd)
   Adios.finalize(rank)
end

-- Run tests
test_1(Mpi.COMM_WORLD)

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
