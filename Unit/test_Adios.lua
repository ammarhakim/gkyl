-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS IO
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
   if rank == 0 then print(msg) end
end

function test_0(comm)
   -- This mimics the test in the adios2 repo:
   --   https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpWriter/helloBPWriter.c
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylad")

   local nx = 10
   local myFloats = Lin.Vec(nx)
   for i = 1, nx do myFloats[i] = 10.*rank + i end

   local shape = {nx*nproc}
   local start = {rank*nx}
   local count = {nx}

   local ad_var = Adios.define_variable(ad_io, "myfloats", Adios.type_double, #shape, shape, start, count, true)

   local ad_engine = Adios.open(ad_io, "myVector.bp", Adios.mode_write)

   local ad_err = Adios.put(ad_engine, ad_var, myFloats:data(), Adios.mode_deferred)

   local _ = Adios.close(ad_engine)

   Adios.finalize(ad)
end

function test_1w(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if nproc > 1 then
      log("Not running test_1 as number of procs > 1")
      return
   end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylad")

   local dim = 2
   local cells = {10, 20}
   local lower = {-math.pi, -math.pi}
   local upper = { math.pi,  math.pi}

   local temp = Lin.Vec(cells[1]) 
   for i = 1, cells[1] do temp[i] = i*1.1 end

   -- Define attributes.
   local cells_attr = Adios.define_attribute_array(ad_io, "numCells", Adios.type_int32_t, cells, dim)
   local lower_attr = Adios.define_attribute_array(ad_io, "lowerBounds", Adios.type_double, lower, dim)
   local upper_attr = Adios.define_attribute_array(ad_io, "upperBounds", Adios.type_double, upper, dim)

   -- Define variables.
   local ad_var_nx   = Adios.define_variable(ad_io, "NX", Adios.type_int32_t, 1, {1}, {0}, {1}, true)
   local ad_var_ny   = Adios.define_variable(ad_io, "NY", Adios.type_int32_t, 1, {1}, {0}, {1}, true)
   local ad_var_temp = Adios.define_variable(ad_io, "temperature", Adios.type_double, 1, {cells[1]}, {0}, {cells[1]}, true)

   local ad_engine = Adios.open(ad_io, "mydata.bp", Adios.mode_write)

   local nxVec, nyVec = Lin.IntVec(1), Lin.IntVec(1)
   nxVec[1], nyVec[1] = cells[1], cells[2]
   local ad_err = Adios.put(ad_engine, ad_var_nx  , nxVec:data(), Adios.mode_deferred)
   local ad_err = Adios.put(ad_engine, ad_var_ny  , nyVec:data(), Adios.mode_deferred)
   local ad_err = Adios.put(ad_engine, ad_var_temp, temp:data(), Adios.mode_deferred)

   local _ = Adios.close(ad_engine)

   Adios.finalize(ad)
end

function test_1rA(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if nproc > 1 then
      log("Not running test_1 as number of procs > 1")
      return
   end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylread")

   local dim = 2
   local cells = {10, 20}
   local lower = {-math.pi, -math.pi}
   local upper = { math.pi,  math.pi}

   -- MF 2023/10/10: Inquiring variables or getting the list of available variables needs
   -- to be done with mode_readRandomAccess or within a Begin/End step block, see issue 3843 in
   -- in the adios2 github: https://github.com/ornladios/ADIOS2/issues/3843
   local ad_engine = Adios.open(ad_io, "mydata.bp", Adios.mode_readRandomAccess)

   local ad_var_nx   = Adios.inquire_variable(ad_io, "NX")
   local ad_var_ny   = Adios.inquire_variable(ad_io, "NY")
   local ad_var_temp = Adios.inquire_variable(ad_io, "temperature")
   local ad_err = Adios.set_selection(ad_var_nx, 1, {0}, {1})
   local ad_err = Adios.set_selection(ad_var_ny, 1, {0}, {1})
   local ad_err = Adios.set_selection(ad_var_nx, 1, {0}, {cells[1]})

   local nxVec, nyVec, temp = Lin.IntVec(1), Lin.IntVec(1), Lin.Vec(cells[1]) 
   local ad_err = Adios.get(ad_engine, ad_var_nx  , nxVec:data(), Adios.mode_sync)
   local ad_err = Adios.get(ad_engine, ad_var_ny  , nyVec:data(), Adios.mode_sync)
   local ad_err = Adios.get(ad_engine, ad_var_temp, temp:data() , Adios.mode_sync)

   local _ = Adios.close(ad_engine)

   assert_equal(10, nxVec[1], "Checking NX")
   assert_equal(20, nyVec[1], "Checking NY")
   for i = 1, cells[1] do
      assert_equal(i*1.1, temp[i], "Checking temperature")
   end
   
   Adios.finalize(ad)
end

function test_1rB(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if nproc > 1 then
      log("Not running test_1 as number of procs > 1")
      return
   end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylread")

   local dim = 2
   local cells = {10, 20}
   local lower = {-math.pi, -math.pi}
   local upper = { math.pi,  math.pi}

   -- MF 2023/10/10: Inquiring variables or getting the list of available variables needs
   -- to be done with mode_readRandomAccess or within a Begin/End step block, see issue 3843 in
   -- in the adios2 github: https://github.com/ornladios/ADIOS2/issues/3843
   local ad_engine = Adios.open(ad_io, "mydata.bp", Adios.mode_read)

   local ad_st_stat = Adios.new_step_status()

   local ad_err = Adios.begin_step(ad_engine, Adios.step_mode_read, -1, ad_st_stat)

   local ad_var_nx   = Adios.inquire_variable(ad_io, "NX")
   local ad_var_ny   = Adios.inquire_variable(ad_io, "NY")
   local ad_var_temp = Adios.inquire_variable(ad_io, "temperature")
   local ad_err = Adios.set_selection(ad_var_nx, 1, {0}, {1})
   local ad_err = Adios.set_selection(ad_var_ny, 1, {0}, {1})
   local ad_err = Adios.set_selection(ad_var_nx, 1, {0}, {cells[1]})

   local nxVec, nyVec, temp = Lin.IntVec(1), Lin.IntVec(1), Lin.Vec(cells[1]) 
   local ad_err = Adios.get(ad_engine, ad_var_nx  , nxVec:data(), Adios.mode_sync)
   local ad_err = Adios.get(ad_engine, ad_var_ny  , nyVec:data(), Adios.mode_sync)
   local ad_err = Adios.get(ad_engine, ad_var_temp, temp:data() , Adios.mode_sync)

   local ad_err = Adios.end_step(ad_engine)

   local _ = Adios.close(ad_engine)

   assert_equal(10, nxVec[1], "Checking NX")
   assert_equal(20, nyVec[1], "Checking NY")
   for i = 1, cells[1] do
      assert_equal(i*1.1, temp[i], "Checking temperature")
   end
   
   Adios.finalize(ad)
end

--function test_2w(comm)
--   local nproc = Mpi.Comm_size(comm)
--   local rank = Mpi.Comm_rank(comm)
--
--   local localSz = 100 -- 100 elements per processor
--   local globalSz = nproc*localSz
--   local myOffset = rank*localSz
--   -- allocate space for local field
--   local temperature = new("double[?]", localSz)
--
--   for i = 0, localSz-1 do
--      temperature[i] = rank*localSz + i
--   end
--
--   Adios.init_noxml(comm)
--   Adios.set_max_buffer_size(16) -- 16 MB chunks
--
--   -- create group and set I/O method
--   local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
--   Adios.select_method(grpId, "MPI", "", "")
--
--   -- define attributes
--   local cells = new("int[1]")
--   cells[0] = globalSz
--   Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, 1, cells)
--
--   local lower = new("double[1]")
--   lower[0] = 0.0;
--   Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, 1, lower)
--
--   local upper = new("double[1]")
--   upper[0] = 1.5;
--   Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, 1, upper)
--
--   -- define variables
--   Adios.define_var(
--      grpId, "temperature", "", Adios.double, tostring(localSz), tostring(globalSz), tostring(myOffset))
--
--   -- open file to write out group
--   local fd = Adios.open("CartField", "adios-test-2.bp", "w", comm)
--   local grpSize = localSz*sizeof("double")
--   local totalSize = Adios.group_size(fd, grpSize)
--   -- write data
--   Adios.write(fd, "temperature", temperature)
--   
--   Adios.close(fd)
--   Adios.finalize(rank)
--end

-- Run tests
test_0(Mpi.COMM_WORLD)
test_1w(Mpi.COMM_WORLD)
test_1rA(Mpi.COMM_WORLD)
test_1rB(Mpi.COMM_WORLD)

--test_2w(Mpi.COMM_WORLD)

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
