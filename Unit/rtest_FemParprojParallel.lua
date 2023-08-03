-- Gkyl ------------------------------------------------------------------------
--
-- Tests for FEM parallel projection operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Lin        = require "Lib.Linalg"
local DecompRegionCalc = require "Lib.CartDecomp"
local Mpi              = require "Comm.Mpi"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local xsys = require "xsys"
local ffi  = require "ffi"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
   "new, copy, fill, sizeof, typeof, metatype")

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then
      print(msg)
   end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

local function createGrid(lo, up, nCells, pDirs, decomp)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
      decomposition = decomp,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   bKind = bKind or "Ser"
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, nGhost)
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {nGhost, nGhost},
      metaData = {polyOrder  = basis:polyOrder(),
                  basisType  = basis:id(),},
   }
   return fld
end

function test_femparproj_1x(nx, polyOrder)
   log("")
   log("Testing 1D parallel FEM projection...")

   local comm = Mpi.COMM_WORLD
   local sz = Mpi.Comm_size(comm)
   if sz < 2 then
      log("Not running test_copyFieldGlobalFromLocal_1D as numProcs less than 2")
      return
   end

   local lower = {-.50}
   local upper = { .50}
   local cells = {nx}
   local numGhost = 1
   log(string.format("nx=%d polyOrder=%d\n", nx, polyOrder))

   local decomp = DecompRegionCalc.CartProd { cuts = {sz}, }
   local grid = createGrid(lower, upper, cells, nil, decomp)

   local noDecomp = DecompRegionCalc.CartProd {
      cuts = { 1 },
      __serTesting = true, -- Hack to create a decomp without cuts.
   }
   local gridGlobal = createGrid(lower, upper, cells, nil, noDecomp)

   local basis = createBasis(grid:ndim(), polyOrder)

   local phi             = createField(grid, basis, numGhost)
   local phiSmooth       = createField(grid, basis, numGhost)
   local phiGlobal       = createField(gridGlobal, basis, numGhost)
   local phiSmoothGlobal = createField(gridGlobal, basis, numGhost)

   local phi_noghost       = createField(grid, basis, 0)
   local phiGlobal_noghost = createField(gridGlobal, basis, 0)

   -- Initialize phi.
   local initPhiFunc = function (t,xn)
      local z = xn[1]
      return math.sin(2.*math.pi*z);
   end
   local initPhi = Updater.ProjectOnBasis {
      onGrid   = grid,  basis = basis,
      evaluate = function (t,xn) return 0. end,  -- Set later.
   }
   initPhi:setFunc(initPhiFunc)
   initPhi:advance(0.,{},{phi})

   -- Copy to a local buffer without ghosts.
   phi_noghost:copyRangeToRange(phi, phi_noghost:localRange(), phi:localRange())

   -- Gather the field from other processes.
   Mpi.Allgather(phi_noghost:dataPointer(), phi_noghost:size(), phi_noghost:elemCommType(),
                 phiGlobal_noghost:dataPointer(), phi_noghost:size(), phiGlobal_noghost:elemCommType(), comm)

   -- Copy to a global field with ghosts
   phiGlobal:copyRangeToRange(phiGlobal_noghost, phiGlobal:localRange(), phiGlobal_noghost:localRange())

   local smoothUpd = Updater.FemParproj {
      onGrid = gridGlobal,  periodicParallelDir = false,
      basis  = basis,       onField             = phiGlobal,
   }
   local t1 = os.clock()
   smoothUpd:advance(0.,{phiGlobal},{phiSmoothGlobal})
   log(string.format("1D parallel smooth took %g s\n", os.clock()-t1))

   -- Copy portion of the global field belonging to this process.
   local lv = phiSmooth:localRange():lowerAsVec()
   local uv = phiSmooth:localRange():upperAsVec()
   local sourceRange = phiSmoothGlobal:globalExtRange():subRange(lv,uv)
   phiSmooth:copyRangeToRange(phiSmoothGlobal, phiSmooth:localRange(), sourceRange)

   local err = createField(grid, basis, numGhost)
   err:combine(1.0, phi, -1.0, phiSmooth)

--   phi:write("phi-rough-1d.bp", 0.0)
--   phiSmooth:write("phi-smooth-1d.bp", 0.0)
--   err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   log(string.format("Average RMS error = %g\n", math.sqrt(lv[1])))
   return math.sqrt(lv[1])
end

function test_femparproj_3x(nx, ny, nz, polyOrder)
   log("")
   log("Testing 3D parallel FEM projection...")

   local comm = Mpi.COMM_WORLD
   local sz = Mpi.Comm_size(comm)
   if sz < 2 then
      log("Not running test_copyFieldGlobalFromLocal_1D as numProcs less than 2")
      return
   end

   local lower = {-2.0, -2.0, -.50}
   local upper = { 2.0,  2.0,  .50}
   local cells = {nx, ny, nz}
   local ncuts = {sz, 1, 1}
   local numGhost = 1
   log(string.format("nx=%d ny=%d nz=%d polyOrder=%d\n", nx, ny, nz, polyOrder))

   local decomp = DecompRegionCalc.CartProd { cuts = ncuts, }
   local grid = createGrid(lower, upper, cells, nil, decomp)

   local noDecomp = DecompRegionCalc.CartProd {
      cuts = { 1, 1, 1 },
      __serTesting = true, -- Hack to create a decomp without cuts.
   }
   local gridGlobal = createGrid(lower, upper, cells, nil, noDecomp)

   local basis = createBasis(grid:ndim(), polyOrder)

   local phi             = createField(grid, basis, numGhost)
   local phiSmooth       = createField(grid, basis, numGhost)
   local phiGlobal       = createField(gridGlobal, basis, numGhost)
   local phiSmoothGlobal = createField(gridGlobal, basis, numGhost)

   local phi_noghost       = createField(grid, basis, 0)
   local phiGlobal_noghost = createField(gridGlobal, basis, 0)

   -- Initialize phi.
   local initPhiFunc = function (t,xn)
      local x, y, z = xn[1], xn[2], xn[3]
      local mu = {0.2, 0.2}
      local sig = 0.3;
      return math.exp(-(math.pow(x-mu[1],2)+math.pow(y-mu[2],2))/(2.0*sig*sig))*math.sin(2.*math.pi*z);
   end
   local initPhi = Updater.ProjectOnBasis {
      onGrid   = grid,  basis = basis,
      evaluate = function (t,xn) return 0. end,  -- Set later.
   }
   initPhi:setFunc(initPhiFunc)
   initPhi:advance(0.,{},{phi})

   -- Copy to a local buffer without ghosts.
   phi_noghost:copyRangeToRange(phi, phi_noghost:localRange(), phi:localRange())

   -- Gather the field from other processes.
   Mpi.Allgather(phi_noghost:dataPointer(), phi_noghost:size(), phi_noghost:elemCommType(),
                 phiGlobal_noghost:dataPointer(), phi_noghost:size(), phiGlobal_noghost:elemCommType(), comm)

   -- Copy to a global field with ghosts
   phiGlobal:copyRangeToRange(phiGlobal_noghost, phiGlobal:localRange(), phiGlobal_noghost:localRange())

--   local smoothUpd = Updater.FemParproj {
--      onGrid = gridGlobal,   periodicParallelDir = false,
--      basis  = basis,        onField             = phiGlobal,
--   }
--   local t1 = os.clock()
--   smoothUpd:advance(0.,{phiGlobal},{phiSmoothGlobal})
--   log(string.format("3D parallel smooth took %g s\n", os.clock()-t1))
   phiSmoothGlobal:copy(phiGlobal)

   -- Copy portion of the global field belonging to this process.
   local lv = phiSmooth:localRange():lowerAsVec()
   local uv = phiSmooth:localRange():upperAsVec()
   local sourceRange = phiSmoothGlobal:globalExtRange():subRange(lv,uv)
   phiSmooth:copyRangeToRange(phiSmoothGlobal, phiSmooth:localRange(), sourceRange)

   local err = createField(grid, basis, numGhost)
   err:combine(1.0, phi, -1.0, phiSmooth)

   phi:write("phi-rough-3d.bp", 0.0)
   phiSmooth:write("phi-smooth-3d.bp", 0.0)
--   err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   log(string.format("Average RMS error = %g\n", math.sqrt(lv[1])))
   return math.sqrt(lv[1])
end

function test_1x_p1()
   log("--- Testing convergence of 1D FEM parallel projection with p=1 ---")
   err1 = test_femparproj_1x(32, 1)
   err2 = test_femparproj_1x(64, 1)
   err3 = test_femparproj_1x(128, 1)
   log(string.format("Order: %g %g\n", err1/err2/2.0, err2/err3/2.0))
   assert_close(4.0, err1/err2/2.0, .05)
   assert_close(4.0, err2/err3/2.0, .05)
   log("")
end

function test_3x_p1()
   log("--- Testing convergence of 3D FEM parallel projection with p=1 ---")
   err1 = test_femparproj_3x(4, 4, 32, 1)
--   err2 = test_femparproj_3x(4, 4, 64, 1)
--   err3 = test_femparproj_3x(4, 4, 128, 1)
--   log(string.format("Order: %g %g\n", err1/err2/2.0, err2/err3/2.0))
--   assert_close(4.0, err1/err2/2.0, .05)
--   assert_close(4.0, err2/err3/2.0, .05)
--   log("")
end

--test_1x_p1()
test_3x_p1()
log("")

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
