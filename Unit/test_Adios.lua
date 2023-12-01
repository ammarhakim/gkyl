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

   local ad_engine = Adios.open(ad_io, "test_Adios_myVector.bp", Adios.mode_write)

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
   local ad_var_nx   = Adios.define_variable(ad_io, "NX", Adios.type_int32_t)
   local ad_var_ny   = Adios.define_variable(ad_io, "NY", Adios.type_int32_t)
   local ad_var_temp = Adios.define_variable(ad_io, "temperature", Adios.type_double, 1, {cells[1]}, {0}, {cells[1]}, true)

   local ad_engine = Adios.open(ad_io, "test_Adios_mydata.bp", Adios.mode_write)

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
   local ad_engine = Adios.open(ad_io, "test_Adios_mydata.bp", Adios.mode_readRandomAccess)

   local ad_var_nx   = Adios.inquire_variable(ad_io, "NX")
   local ad_var_ny   = Adios.inquire_variable(ad_io, "NY")
   local ad_var_temp = Adios.inquire_variable(ad_io, "temperature")
   local ad_err = Adios.set_selection(ad_var_temp, 1, {0}, {cells[1]})

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
   local ad_engine = Adios.open(ad_io, "test_Adios_mydata.bp", Adios.mode_read)

   local ad_st_stat = Adios.new_step_status()

   local ad_err = Adios.begin_step(ad_engine, Adios.step_mode_read, -1, ad_st_stat)

   local ad_var_nx   = Adios.inquire_variable(ad_io, "NX")
   local ad_var_ny   = Adios.inquire_variable(ad_io, "NY")
   local ad_var_temp = Adios.inquire_variable(ad_io, "temperature")
   local ad_err = Adios.set_selection(ad_var_temp, 1, {0}, {cells[1]})

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

function test_2w(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   local localSz = 10 -- 10 elements per processor
   -- Allocate and initialize local field.
   local globalSz, myOffset = nproc*localSz, rank*localSz
   local temp = Lin.Vec(localSz)
   for i = 1, localSz do temp[i] = rank*localSz + i end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylad")

   local dim = 1
   local cells = {globalSz}
   local lower = {-math.pi}
   local upper = { math.pi}

   -- Define attributes.
   local cells_attr = Adios.define_attribute_array(ad_io, "numCells", Adios.type_int32_t, cells, dim)
   local lower_attr = Adios.define_attribute_array(ad_io, "lowerBounds", Adios.type_double, lower, dim)
   local upper_attr = Adios.define_attribute_array(ad_io, "upperBounds", Adios.type_double, upper, dim)

   -- Define variables.
   local ad_var_temp = Adios.define_variable(ad_io, "temperature", Adios.type_double, 1, {globalSz}, {myOffset}, {localSz}, true)

   local ad_engine = Adios.open(ad_io, "test_Adios_mydata.bp", Adios.mode_write)

   local ad_err = Adios.put(ad_engine, ad_var_temp, temp:data(), Adios.mode_deferred)

   local _ = Adios.close(ad_engine)

   Adios.finalize(ad)
end

function test_2rA(comm)
   local nproc, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   local localSz = 10 -- 10 elements per processor
   -- Allocate and initialize local field.
   local globalSz, myOffset = nproc*localSz, rank*localSz
   local mom = {temperature = Lin.Vec(localSz)}
   for i = 1, localSz do mom.temperature[i] = rank*localSz + i end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylread")

   local dim = 1
   local grid = { numCells    = {globalSz},
                  lowerBounds = {-math.pi},
                  upperBounds = { math.pi}, }

   -- MF 2023/10/10: Inquiring variables or getting the list of available variables needs
   -- to be done with mode_readRandomAccess or within a Begin/End step block, see issue 3843 in
   -- in the adios2 github: https://github.com/ornladios/ADIOS2/issues/3843
   local ad_engine = Adios.open(ad_io, "test_Adios_mydata.bp", Adios.mode_readRandomAccess)

   local attr_names = Adios.available_attributes(ad_io)
   local ad_attrs, attrs = {}, {}
   for i, nm in ipairs(attr_names) do
      ad_attrs[i] = Adios.inquire_attribute(ad_io, nm)
      local name = Adios.attribute_name(ad_attrs[i])
      assert_equal(nm, name, "test_2rA: checking attribute names")

      attrs[nm], sz = Adios.attribute_data(ad_attrs[i])

      assert_equal(grid[nm][1], attrs[nm][1], "test_2rA: checking attribute values")
   end

   local var_names = Adios.available_variables(ad_io)
   local ad_vars, vars = {}, {}
   for i, nm in ipairs(var_names) do
      ad_vars[i] = Adios.inquire_variable(ad_io, nm)
      local name = Adios.variable_name(ad_vars[i])
      assert_equal(nm, name, "test_2rA: checking variable names")

      local ndims = Adios.variable_ndims(ad_vars[i])
      assert_equal(ndims, 1, "test_2rA: checking variable ndims")

      local shape = Adios.variable_shape(ad_vars[i])
      assert_equal(shape[1], globalSz, "test_2rA: checking variable shape")

      local ad_err = Adios.set_selection(ad_vars[i], 1, {rank*localSz}, {localSz})
      vars[nm] = Lin.Vec(localSz)
      local ad_err = Adios.get(ad_engine, ad_vars[i], vars[nm]:data() , Adios.mode_sync)
   end

   local _ = Adios.close(ad_engine)

   for i = 1, localSz do
      assert_equal(mom.temperature[i], vars.temperature[i], "test_2rA: checking variable values")
   end
   
   Adios.finalize(ad)
end

function test_3w(comm)
   -- Test writing using 2 different communicators.
   -- Simulate 2 species, distributing a (velocity) moment across 2 configuration space ranks and 2 species ranks.
   local defaultCommSize = 4
   local worldSize, worldRank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if worldSize ~= defaultCommSize then
      log(string.format("Not running test_3w as comm size not exactly %d", defaultCommSize))
      return
   end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylio")

   local numCutsConf  = 2

   local confColor    = math.floor(worldRank / numCutsConf)
   local confComm     = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, worldRank)
   local speciesColor = worldRank % numCutsConf
   local speciesComm  = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, worldRank)
   local confSize, confRank       = Mpi.Comm_size(confComm),    Mpi.Comm_rank(confComm)
   local speciesSize, speciesRank = Mpi.Comm_size(speciesComm), Mpi.Comm_rank(speciesComm)

   local dim = 1
   local cells = {12}
   local lower = {-math.pi}
   local upper = { math.pi}

   -- Allocate data distributed across different comms.
   local cellsLocal = cells[1]/confSize
   local moms = {}
   for s = 1, speciesSize do
      moms[s] = Lin.Vec(cellsLocal)
      for i = 1, cellsLocal do moms[s][i] = 1.1*(confRank*cellsLocal+i-1)+(s-1)*10. end
   end

   -- These definitions can also go inside of the for loop and if-statement below.
   -- Define attributes.
   local cells_attr = Adios.define_attribute_array(ad_io, "numCells", Adios.type_int32_t, cells, dim)
   local lower_attr = Adios.define_attribute_array(ad_io, "lowerBounds", Adios.type_double, lower, dim)
   local upper_attr = Adios.define_attribute_array(ad_io, "upperBounds", Adios.type_double, upper, dim)
   -- Define variables.
   local ad_var_mom = Adios.define_variable(ad_io, "mom", Adios.type_double, 1, cells, {confRank*cellsLocal}, {cellsLocal}, true)

   for s = 1, speciesSize do
      if speciesRank == s-1 then
         local ad_engine = Adios.open_new_comm(ad_io, "test_Adios_mom_"..tostring(s-1)..".bp", Adios.mode_write, confComm)
   
         local ad_err = Adios.put(ad_engine, ad_var_mom, moms[s]:data(), Adios.mode_deferred)
   
         local _ = Adios.close(ad_engine)
      end
   end

   -- Test the case in which only one species rank writes a file.
   field = Lin.Vec(cellsLocal)
   for i = 1, cellsLocal do field[i] = 3.*(confRank*cellsLocal+i-1) end

   if speciesRank == 0 then
      -- Define variables.
      local ad_var_field = Adios.define_variable(ad_io, "field", Adios.type_double, 1, cells, {confRank*cellsLocal}, {cellsLocal}, true)
      local ad_engine = Adios.open_new_comm(ad_io, "test_Adios_field.bp", Adios.mode_write, confComm)
      local ad_err = Adios.put(ad_engine, ad_var_field, field:data(), Adios.mode_deferred)
      local _ = Adios.close(ad_engine)
   end

   Adios.finalize(ad)
end

function test_3r(comm)
   -- Test reading using 2 different communicators.
   -- Simulate 2 species, distributing a (velocity) moment across 2 configuration space ranks and 2 species ranks.
   local defaultCommSize = 4
   local worldSize, worldRank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if worldSize ~= defaultCommSize then
      log(string.format("Not running test_3r as comm size not exactly %d", defaultCommSize))
      return
   end

   local ad = Adios.init_mpi(comm)
   local ad_io = Adios.declare_io(ad, "gkylio")

   local numCutsConf  = 2

   local confColor    = math.floor(worldRank / numCutsConf)
   local confComm     = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, worldRank)
   local speciesColor = worldRank % numCutsConf
   local speciesComm  = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, worldRank)
   local confSize, confRank       = Mpi.Comm_size(confComm),    Mpi.Comm_rank(confComm)
   local speciesSize, speciesRank = Mpi.Comm_size(speciesComm), Mpi.Comm_rank(speciesComm)

   local dim = 1
   local cells = {12}
   local lower = {-math.pi}
   local upper = { math.pi}

   -- Allocate data distributed across different comms.
   local cellsLocal = cells[1]/confSize
   local moms, momsRef = {}, {}
   for s = 1, speciesSize do
      moms[s], momsRef[s] = Lin.Vec(cellsLocal), Lin.Vec(cellsLocal)
      for i = 1, cellsLocal do momsRef[s][i] = 1.1*(confRank*cellsLocal+i-1)+(s-1)*10. end
   end

   for s = 1, speciesSize do
      if speciesRank == s-1 then
         local ad_engine = Adios.open_new_comm(ad_io, "test_Adios_mom_"..tostring(s-1)..".bp", Adios.mode_readRandomAccess, confComm)
         local ad_var_mom = Adios.inquire_variable(ad_io, "mom")
         local ad_err = Adios.set_selection(ad_var_mom, 1, {confRank*cellsLocal}, {cellsLocal})
   
         local ad_err = Adios.get(ad_engine, ad_var_mom, moms[s]:data(), Adios.mode_deferred)
   
         local _ = Adios.close(ad_engine)

         -- Check data read.
         for i = 1, cellsLocal do 
            assert_equal(momsRef[s][i], moms[s][i], "test_3r: checking mom data read.")
         end
      end
   end

   -- Read the field from all species ranks.
   field, fieldRef = Lin.Vec(cellsLocal), Lin.Vec(cellsLocal)
   for i = 1, cellsLocal do fieldRef[i] = 3.*(confRank*cellsLocal+i-1) end

   local ad_engine = Adios.open_new_comm(ad_io, "test_Adios_field.bp", Adios.mode_readRandomAccess, confComm)
   local ad_var_field = Adios.inquire_variable(ad_io, "field")
   local ad_err = Adios.set_selection(ad_var_field, 1, {confRank*cellsLocal}, {cellsLocal})
   local ad_err = Adios.get(ad_engine, ad_var_field, field:data(), Adios.mode_deferred)
   local _ = Adios.close(ad_engine)

   -- Check data read.
   for i = 1, cellsLocal do 
      assert_equal(fieldRef[i], field[i], "test_3r: checking field data read.")
   end


   Adios.finalize(ad)
end

-- Run tests
test_0(Mpi.COMM_WORLD)
test_1w(Mpi.COMM_WORLD)
test_1rA(Mpi.COMM_WORLD)
test_1rB(Mpi.COMM_WORLD)
test_2w(Mpi.COMM_WORLD)
test_2rA(Mpi.COMM_WORLD)
test_3w(Mpi.COMM_WORLD)
test_3r(Mpi.COMM_WORLD)

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
