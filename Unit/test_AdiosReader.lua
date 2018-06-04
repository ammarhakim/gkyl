-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS IO
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosReader = require "Io.AdiosReader"
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

function test_1(comm)
   local reader = AdiosReader.Reader("t2-two-stream_elc_10.bp")

   -- check if expected variables are present in file
   assert_equal(true, reader:hasVar("time"), "Checking for 'time'")
   assert_equal(true, reader:hasVar("frame"), "Checking for 'frame'")
   assert_equal(true, reader:hasVar("CartGridField"), "Checking for 'CartGridField'")
   assert_equal(false, reader:hasVar("Nothing"), "Checking for 'Nothing'")

   local time = reader:getVar("time")
   assert_equal("time", time.name, "Checking time ")
   assert_equal(50.0, time:read(), "Checking time-value")
   assert_equal("double", time.type, "Checking time-type")   

   local frame = reader:getVar("frame")
   assert_equal("frame", frame.name, "Checking frame ")
   assert_equal(10, frame:read(), "Checking frame-value")
   assert_equal("integer", frame.type, "Checking frame-type")   

   local field = reader:getVar("CartGridField")
   assert_equal("CartGridField", field.name, "Checking name")
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

-- Run tests
test_1(Mpi.COMM_WORLD)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
