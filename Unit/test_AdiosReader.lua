-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS IO
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosReader = require "Io.AdiosReader"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Range = require "Range"
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
   local reader = AdiosReader.Reader("../Unit/t2-two-stream_elc_10.bp", comm)

   -- check if expected variables are present in file
   assert_equal(true, reader:hasVar("time"), "Checking for 'time'")
   assert_equal(true, reader:hasVar("frame"), "Checking for 'frame'")
   assert_equal(true, reader:hasVar("CartGridField"), "Checking for 'CartGridField'")
   assert_equal(false, reader:hasVar("Nothing"), "Checking for 'Nothing'")

   local time = reader:getVar("time")
   assert_equal("time", time.name, "Checking time ")
   assert_equal(50.001874949619, time:read(), "Checking time-value")
   assert_equal("double", time.type, "Checking time-type")

   local frame = reader:getVar("frame")
   assert_equal("frame", frame.name, "Checking frame ")
   assert_equal(10, frame:read(), "Checking frame-value")
   assert_equal("integer", frame.type, "Checking frame-type")

   local field = reader:getVar("CartGridField")
   assert_equal("CartGridField", field.name, "Checking name")
   local data = field:read()

   -- check if expected attributes are present in file
   assert_equal(true, reader:hasAttr("numCells"), "Num cells")
   assert_equal(true, reader:hasAttr("lowerBounds"), "lower bounds")
   assert_equal(true, reader:hasAttr("upperBounds"), "upper bounds")
   assert_equal(true, reader:hasAttr("grid"), "grid")
   assert_equal(true, reader:hasAttr("type"), "type")
   assert_equal(false, reader:hasAttr("Nothing"), "checking 'Nothing'")

   -- read attributes
   local numCells = reader:getAttr("numCells"):read()
   assert_equal(64, numCells[1], "Number of cells in X")
   assert_equal(32, numCells[2], "Number of cells in y")

   local lowerBounds = reader:getAttr("lowerBounds"):read()
   assert_equal(-2*math.pi, lowerBounds[1], "Lower bound in X")
   assert_equal(-6, lowerBounds[2], "Lower bound in y")

   local upperBounds = reader:getAttr("upperBounds"):read()
   assert_equal(2*math.pi, upperBounds[1], "Upper bound in X")
   assert_equal(6, upperBounds[2], "Upper bound in y")

   local gtype = reader:getAttr("type"):read()
   assert_equal("uniform", gtype[1], "Grid type")

   local gfile = reader:getAttr("grid"):read()
   assert_equal("grid", gfile[1], "Grid file name")

   reader:close()
end

function test_2(comm)
   local reader = AdiosReader.Reader("../Unit/test_DynVector_test_1.bp", comm)

   assert_equal(true, reader:hasVar("time"), "Checking for 'time'")
   assert_equal(true, reader:hasVar("frame"), "Checking for 'frame'")

   local timeMesh = reader:getVar("TimeMesh"):read()
   local dataVar = reader:getVar("Data")
   local data = dataVar:read()

   local dynRange = Range.Range({1, 1}, {dataVar.shape[1], dataVar.shape[2]})
   local indexer = Range.makeRowMajorIndexer(dynRange)

   for i = 1, dynRange:shape(1) do
      assert_equal(2.5*i^2,  data[indexer(i,1)], "Testing read data")
      assert_equal(2.5*i^2+0.5,  data[indexer(i,2)], "Testing read data")
   end

   reader:close()
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

-- Run tests
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
