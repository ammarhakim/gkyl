-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute integrated quantities.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Messenger  = require "Comm.Messenger"
local DecompRegionCalc = require "Lib.CartDecomp"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

local function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then print(msg) end
end

local function test_ser_op_none()
   local comm = Mpi.COMM_WORLD
   local sz, rnk = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 1 then
      log("test_ser_op_none: Not running as number of procs not exactly 1.")
      return
   end

   local pOrder = 2

   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = pOrder }
   local confBasis  = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = pOrder }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid = phaseGrid,  numComponents = phaseBasis:numBasis(),
      ghost  = {0, 0},
   }
   local numDensity = DataStruct.Field {
      onGrid = confGrid,  numComponents = confBasis:numBasis(),
      ghost  = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,  basis = phaseBasis,
      evaluate = function (t, xn)
         return 1/((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},  phaseBasis = phaseBasis,
      onGrid      = phaseGrid,                confBasis  = confBasis,
      moment      = "M0",
   }
   calcNumDensity:advance(0.0, {distf}, {numDensity})

   -- Compute integrated f and integrated moment.
   local intF = DataStruct.DynVector { numComponents = 1, }
   local intM0 = DataStruct.DynVector { numComponents = 1, }
   local intQuantP = Updater.CartFieldIntegrate {
      onGrid = phaseGrid,   numComponents = 1,
      basis  = phaseBasis,
   }
   local intQuantC = Updater.CartFieldIntegrate {
      onGrid = confGrid,   numComponents = 1,
      basis  = confBasis,
   }

   intQuantP:advance(0.0, {distf}, {intF})
   intQuantC:advance(0.0, {numDensity}, {intM0})

   local _, NP = intF:lastData()
   local _, NC = intM0:lastData()

   assert_equal(1, NP[1], "test_ser_op_none: Checking phase moment.")
   assert_equal(1, NC[1], "test_ser_op_none: Checking conf moment.")
end

local function test_ser_op_sq()
   local comm = Mpi.COMM_WORLD
   local sz, rnk = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 1 then
      log("test_ser_op_sq: Not running as number of procs not exactly 1.")
      return
   end

   local pOrder = 2

   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = pOrder }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid = phaseGrid,  numComponents = phaseBasis:numBasis(),
      ghost  = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,  basis = phaseBasis,
      evaluate = function (t, xn)
         return 1/math.sqrt((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Compute integrated f and integrated moment.
   local fL2 = DataStruct.DynVector { numComponents = 1, }
   local intQuantP = Updater.CartFieldIntegrate {
      onGrid = phaseGrid,   numComponents = 1,
      basis  = phaseBasis,  operator      = 'sq',
   }

   intQuantP:advance(0.0, {distf}, {fL2})

   local _, intL2 = fL2:lastData()

   assert_equal(1, intL2[1], "test_ser_op_sq: Checking moment.")
end

local function test_ser_op_abs()
   local comm = Mpi.COMM_WORLD
   local sz, rnk = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 1 then
      log("test_ser_op_abs: Not running as number of procs not exactly 1.")
      return
   end

   local pOrder = 2

   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = pOrder }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid = phaseGrid,  numComponents = phaseBasis:numBasis(),
      ghost  = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,  basis = phaseBasis,
      evaluate = function (t, xn)
         return -1/((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Compute integrated f and integrated moment.
   local absF = DataStruct.DynVector { numComponents = 1, }
   local intQuantP = Updater.CartFieldIntegrate {
      onGrid = phaseGrid,   numComponents = 1,
      basis  = phaseBasis,  operator      = 'abs',
   }

   intQuantP:advance(0.0, {distf}, {absF})

   local _, intAbsF = absF:lastData()

   assert_equal(1, intAbsF[1], "test_ser_op_abs: Checking moment.")
end

local function test_ser_op_none_par()
   local comm = Mpi.COMM_WORLD
   local sz, rnk = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("test_ser_op_none_par: Not running as number of procs not exactly 2.")
      return
   end

   local pOrder = 1
   local phaseLower = {0., -6.0, -3.0}
   local phaseUpper = {1.,  6.0,  3.0}
   local phaseCells = {4, 16, 12}
   local decompCuts = {2}

   local commManager = Messenger{
      cells      = {phaseCells[1]},   decompCutsConf = decompCuts,
      numSpecies = 1,  parallelizeSpecies = false,
   }
   -- Phase-space and config-space grids.
   local confGrid = Grid.RectCart {
      lower = { phaseLower[1], },  decomposition = commManager:getConfDecomp(),
      upper = { phaseUpper[1], },  messenger     = commManager,
      cells = { phaseCells[1] },
   }
   local decompCutsPhase = {decompCuts[1], 1, 1} 
   local decompPhase = DecompRegionCalc.CartProd {
      cuts = decompCutsPhase,  comm = confGrid:commSet().comm,
   }
   local phaseGrid = Grid.RectCart {
      lower = phaseLower,  decomposition = decompPhase,
      upper = phaseUpper,  messenger     = confGrid:getMessenger(),
      cells = phaseCells,
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = pOrder }
   local confBasis  = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = pOrder }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid = phaseGrid,  numComponents = phaseBasis:numBasis(),
      ghost  = {1, 1},
   }
   local numDensity = DataStruct.Field {
      onGrid = confGrid,  numComponents = confBasis:numBasis(),
      ghost  = {1, 1},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,  basis = phaseBasis,
      evaluate = function (t, xn)
         local f = 1.
         for d = 1,phaseBasis:ndim() do f = f*1/(phaseGrid:upper(d)-phaseGrid:lower(d)) end
         return f
      end
   }
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},  phaseBasis = phaseBasis,
      onGrid      = phaseGrid,                confBasis  = confBasis,
      moment      = "M0",
   }
   calcNumDensity:advance(0.0, {distf}, {numDensity})

   -- Compute integrated f and integrated moment.
   local intF = DataStruct.DynVector { numComponents = 1, }
   local intM0 = DataStruct.DynVector { numComponents = 1, }
   local intQuantP = Updater.CartFieldIntegrate {
      onGrid = phaseGrid,   numComponents = 1,
      basis  = phaseBasis,
   }
   local intQuantC = Updater.CartFieldIntegrate {
      onGrid = confGrid,   numComponents = 1,
      basis  = confBasis,
   }

   intQuantP:advance(0.0, {distf}, {intF})
   intQuantC:advance(0.0, {numDensity}, {intM0})

   local _, NP = intF:lastData()
   local _, NC = intM0:lastData()

   assert_equal(1, NP[1], "test_ser_op_none: Checking phase moment.")
   assert_equal(1, NC[1], "test_ser_op_none: Checking conf moment.")
end

test_ser_op_none()
test_ser_op_sq()
test_ser_op_abs()
test_ser_op_none_par()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
