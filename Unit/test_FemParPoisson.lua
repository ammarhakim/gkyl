-- Gkyl ------------------------------------------------------------------------
--
-- Tests for FEM parallel Poisson solver
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Basis = require "Basis"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local CalcDiagnostic = require "Updater.CalcDiagnostic"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

function test_solve1d(nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 1D parallel Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     bcBack = { T = "D", V = 0.0 },
     bcFront = { T = "D", V = 0.0 },
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("1D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, -12
   local c1 = 0
   local c0 = -a/12 - 1/2 - b/30
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[1]
                    return  1+a*z^2+b*z^4
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[1]
                    return z^2/2+a*z^4/12+b*z^6/30+c0*z+c1
                 end
   }
   initExactSolModal:advance(0.,0.,{},{exactSolModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   print("Solving...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("1D parallel Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("phi-solution-1d.bp", 0.0)
   --exactSolModal:write("exact-solution-1d.bp", 0.0)
   --err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_smooth1d(nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 1D parallel Poisson smooth...")
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     --bcBack = { T = "N", V = 0.0 },
     --bcFront = { T = "N", V = 0.0 },
     periodicDirs = {1},
     laplacianWeight = 0.0,
     modifierConstant = 1.0,
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("1D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[1]
                    return math.cos(2*math.pi*2*z)
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   print("Smoothing...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("1D parallel Poisson smooth took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, srcModal, -1.0, phiModal)

   --phiModal:write("phi-solution-1d.bp", 0.0)
   --exactSolModal:write("exact-solution-1d.bp", 0.0)
   --err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_solve2d(nx, nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 2D parallel Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {nx, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   io.write("nx=",nx," nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     bcBack = { T = "D", V = 0.0 },
     bcFront = { T = "D", V = 0.0 },
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("2D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, -12
   local c1 = 0
   local c0 = -a/12 - 1/2 - b/30
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local z = xn[2]
                    return  (1+a*z^2+b*z^4)*(1)
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local z = xn[2]
                    return (z^2/2+a*z^4/12+b*z^6/30+c0*z+c1)*(1)
                 end
   }
   initExactSolModal:advance(0.,0.,{},{exactSolModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   print("Solving...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("2D parallel Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("phi-solution-2d.bp", 0.0)
   --exactSolModal:write("exact-solution-2d.bp", 0.0)
   --err:write("error-2d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_smooth2d(nx, nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 2D parallel Poisson smooth...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {nx, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   io.write("nx=",nx," nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     bcBack = { T = "N", V = 0.0 },
     bcFront = { T = "N", V = 0.0 },
     laplacianWeight = 0.0,
     modifierConstant = 1.0,
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("2D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, -12
   local c1 = 0
   local c0 = -a/12 - 1/2 - b/30
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local z = xn[2]
                    return  (1+a*z^2+b*z^4)*(1)
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   print("Smoothing...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("2D parallel Poisson smooth took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, srcModal, -1.0, phiModal)

   --phiModal:write("phi-solution-2d.bp", 0.0)
   --exactSolModal:write("exact-solution-2d.bp", 0.0)
   --err:write("error-2d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_solve3d(nx, ny, nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 3D parallel Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {1.0, 1.0, 1.0},
      cells = {nx, ny, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 3, polyOrder = p }
   io.write("nx=",nx," ny=",ny," nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     bcBack = { T = "D", V = 0.0 },
     bcFront = { T = "D", V = 0.0 },
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("3D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, -12
   local c1 = 0
   local c0 = -a/12 - 1/2 - b/30
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[3]
                    return  1+a*z^2+b*z^4
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[3]
                    return z^2/2+a*z^4/12+b*z^6/30+c0*z+c1
                 end
   }
   initExactSolModal:advance(0.,0.,{},{exactSolModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   print("Solving...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("3D parallel Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("phi-solution-1d.bp", 0.0)
   --exactSolModal:write("exact-solution-1d.bp", 0.0)
   --err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_smooth3d(nx, ny, nz, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 3D parallel Poisson smooth...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {1.0, 1.0, 1.0},
      cells = {nx, ny, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 3, polyOrder = p }
   io.write("nx=",nx," ny=",ny," nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParPoisson {
     onGrid = grid,
     basis = basis,
     bcBack = { T = "N", V = 0.0 },
     bcFront = { T = "N", V = 0.0 },
     laplacianWeight = 0.0,
     modifierConstant = 1.0,
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("3D parallel Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, -12
   local c1 = 0
   local c0 = -a/12 - 1/2 - b/30
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local z = xn[3]
                    return  1+a*z^2+b*z^4
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   phiModal:copy(srcModal)
   print("Smoothing...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{phiModal},{phiModal})
   local t2 = os.clock()
   io.write("3D parallel Poisson smooth took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, srcModal, -1.0, phiModal)

   --phiModal:write("phi-solution-1d.bp", 0.0)
   --exactSolModal:write("exact-solution-1d.bp", 0.0)
   --err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_solve1d_p1()
  print("--- Testing convergence of 1D solver with p=1 ---")
  err1 = test_solve1d(32, 1)
  err2 = test_solve1d(64, 1)
  err3 = test_solve1d(128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(2.0, err1/err2/2.0, .01)
  assert_close(2.0, err2/err3/2.0, .01)
  print()
end

function test_solve1d_p2()
  print("--- Testing convergence of 1D solver with p=2 ---")
  err1 = test_solve1d(32, 2)
  err2 = test_solve1d(64, 2)
  err3 = test_solve1d(128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .01)
  assert_close(4.0, err2/err3/2.0, .01)
  print()
end

function test_solve2d_p1()
  print("--- Testing convergence of 2D solver with p=1 ---")
  err1 = test_solve2d(1, 32, 1)
  err2 = test_solve2d(1, 64, 1)
  err3 = test_solve2d(1, 128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(2.0, err1/err2/2.0, .01)
  assert_close(2.0, err2/err3/2.0, .01)
  print()
end

function test_solve2d_p2()
  print("--- Testing convergence of 2D solver with p=2 ---")
  err1 = test_solve2d(4, 32, 2)
  err2 = test_solve2d(4, 64, 2)
  err3 = test_solve2d(4, 128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .01)
  assert_close(4.0, err2/err3/2.0, .01)
  print()
end

function test_solve3d_p1()
  print("--- Testing convergence of 3D solver with p=1 ---")
  err1 = test_solve3d(4, 4, 32, 1)
  err2 = test_solve3d(4, 4, 64, 1)
  err3 = test_solve3d(4, 4, 128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(2.0, err1/err2/2.0, .01)
  assert_close(2.0, err2/err3/2.0, .01)
  print()
end

function test_solve3d_p2()
  print("--- Testing convergence of 3D solver with p=2 ---")
  err1 = test_solve3d(4, 4, 32, 2)
  err2 = test_solve3d(4, 4, 64, 2)
  err3 = test_solve3d(4, 4, 128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .01)
  assert_close(4.0, err2/err3/2.0, .01)
  print()
end

function test_smooth1d_p1()
  print("--- Testing convergence of 1D smoother with p=1 ---")
  err1 = test_smooth1d(32, 1)
  err2 = test_smooth1d(64, 1)
  err3 = test_smooth1d(128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .05)
  assert_close(4.0, err2/err3/2.0, .05)
  print()
end

function test_smooth1d_p2()
  print("--- Testing convergence of 1D smoother with p=2 ---")
  err1 = test_smooth1d(32, 2)
  err2 = test_smooth1d(64, 2)
  err3 = test_smooth1d(128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .1)
  assert_close(4.0, err2/err3/2.0, .1)
  print()
end

function test_smooth2d_p1()
  print("--- Testing convergence of 2D smoother with p=1 ---")
  err1 = test_smooth2d(1, 32, 1)
  err2 = test_smooth2d(1, 64, 1)
  err3 = test_smooth2d(1, 128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .05)
  assert_close(4.0, err2/err3/2.0, .05)
  print()
end

function test_smooth2d_p2()
  print("--- Testing convergence of 2D smoother with p=2 ---")
  err1 = test_smooth2d(4, 32, 2)
  err2 = test_smooth2d(4, 64, 2)
  err3 = test_smooth2d(4, 128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .1)
  assert_close(4.0, err2/err3/2.0, .1)
  print()
end

function test_smooth3d_p1()
  print("--- Testing convergence of 3D smoother with p=1 ---")
  err1 = test_smooth3d(4, 4, 32, 1)
  err2 = test_smooth3d(4, 4, 64, 1)
  err3 = test_smooth3d(4, 4, 128, 1)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .05)
  assert_close(4.0, err2/err3/2.0, .05)
  print()
end

function test_smooth3d_p2()
  print("--- Testing convergence of 3D smoother with p=2 ---")
  err1 = test_smooth3d(4, 4, 32, 2)
  err2 = test_smooth3d(4, 4, 64, 2)
  err3 = test_smooth3d(4, 4, 128, 2)
  print("Order:", err1/err2/2.0, err2/err3/2.0)
  assert_close(4.0, err1/err2/2.0, .1)
  assert_close(4.0, err2/err3/2.0, .1)
  print()
end

local t1 = os.clock()
test_solve1d_p1()
test_solve1d_p2()
test_solve2d_p1()
test_solve2d_p2()
test_solve3d_p1()
test_solve3d_p2()
test_smooth1d_p1()
test_smooth1d_p2()
test_smooth2d_p1()
test_smooth2d_p2()
test_smooth3d_p1()
test_smooth3d_p2()
local t2 = os.clock()
print()
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
io.write("Total test time: ", t2-t1, " s\n")
