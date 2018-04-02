-- Gkyl ------------------------------------------------------------------------
--
-- Tests for 2D FEM Poisson solver
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

function test_init_nonperiodic()
   local nx = 256
   local ny = 256
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {nx, ny},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local t1 = os.clock()
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     bcLeft = { T = "D", V = 0.0 },
     bcRight = { T = "D", V = 0.0 },
     bcBottom = { T = "N", V = 0.0 },
     bcTop = { T = "D", V = 0.0 },
     -- periodicDirs = {1,2},
     --writeStiffnessMatrix = true,
   }
   local t2 = os.clock()
   print("256x256 Poisson init took", t2-t1, "s")
   if nx==256 and ny==256 and t2-t1>5 then
     print("Warning... Poisson solver initialization took longer than it should...")
   end

   assert_equal(nx, poisson._nx)
   assert_equal(ny, poisson._ny)
   assert_equal(1.0/256, poisson._dx, "dx")
   assert_equal(0, poisson:bcType(0,0))
   assert_equal(0, poisson:bcType(0,1))
   assert_equal(1, poisson:bcType(1,0))
   assert_equal(0, poisson:bcType(1,1))
   assert_equal(0.0, poisson:bcValue(0,0))
   assert_equal(0.0, poisson:bcValue(0,1))
   assert_equal(0.0, poisson:bcValue(1,0))
   assert_equal(0.0, poisson:bcValue(1,1))
end

function test_init_allperiodic()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {4, 8},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1,2},
   }
   assert_equal(true, poisson._allPeriodic, "allPeriodic")
   assert_equal(true, poisson._isDirPeriodic[0], "periodic dir 0")
   assert_equal(true, poisson._isDirPeriodic[1], "periodic dir 1")
end

function test_init_someperiodic()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {4, 8},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1},
     bcBottom = { T = "N", V = 0.0 },
     bcTop = { T = "D", V = 0.0 },
   }
   assert_equal(false, poisson._allPeriodic, "allPeriodic")
   assert_equal(true,  poisson._isDirPeriodic[0], "periodic dir 0")
   assert_equal(false, poisson._isDirPeriodic[1], "periodic dir 1")
end

function test_init_errorcheck()
   assert_equal(false, pcall(Updater.FemPerpPoisson, {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1},
     bcBottom = { T = "N", V = 0.0 },
     bcLeft = { T = "D", V = 0.0 },
   }), "mismatched BCs")

   assert_equal(false, pcall(Updater.FemPerpPoisson, {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1},
     bcBottom = { T = "N", V = 0.0 },
     bcTop = { T = "V", V = 0.0 },
   }), "bad BC value")

   assert_equal(false, pcall(Updater.FemPerpPoisson, {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1},
     bcBottom = { T = "N", V = 0.0 },
   }), "not enough BCs")
end

function test_solve2d(nx, ny, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 2D Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {nx, ny},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   io.write("nx=",nx," ny=", ny," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     bcLeft = { T = "D", V = 0.0 },
     bcRight = { T = "D", V = 0.0 },
     bcBottom = { T = "N", V = 0.0 },
     bcTop = { T = "D", V = 0.0 },
     -- periodicDirs = {1,2},
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("2D Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local a, b = 2, 5
   local c1, d0 = 0, 0
   local c0 = a/12 - 1/2
   local d1 = b/12 - 1/2
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local y = xn[2]
                    local t1 = (1-a*x^2)*(-b*y^4/12 + y^2/2 + d0*y + d1)
                    local t2 = (1-b*y^2)*(-a*x^4/12 + x^2/2 + c0*x + c1)
                    return t1+t2
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
                    local y = xn[2]
                    local t1 = x^2/2-a*x^4/12+c0*x+c1
                    local t2 = y^2/2-b*y^4/12+d0*y+d1
                    return t1*t2
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
   io.write("2D Poisson solve took total of ", t2-t1, " s\n")

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

function test_solve2d_periodic(nx, ny, p)
   print()
   print("Testing 2D periodic Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {2*math.pi, 2*math.pi},
      cells = {nx, ny},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   io.write("nx=",nx," ny=", ny," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1,2},
     writeStiffnessMatrix = false,
   }
   local t2 = os.clock()
   io.write("2D periodic Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   -- initialization constants
   local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
   local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                 local x = xn[1]
                 local y = xn[2]
                 local t1, t2 = 0.0, 0.0
                 local f = 0.0
                 for m = 0,2 do
                    for n = 0,2 do
                       t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
                       t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
                       f = f + -(m*m+n*n)*(t1+t2)
                    end
                 end
                 return f/50.0
              end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
   }
   local initfunc = function (t,xn)
              local x = xn[1]
              local y = xn[2]
              local t1, t2 = 0.0, 0.0
              local f = 0.0
              for m = 0,2 do
                 for n = 0,2 do
                    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
                    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
                    f = f + (t1+t2)
                 end
              end
              f = f/50.0
              return f-.6
           end
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = initfunc,
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
   io.write("2D periodic Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("periodic-phi-solution-2d.bp", 0.0)
   --exactSolModal:write("periodic-exact-solution-2d.bp", 0.0)
   --err:write("periodic-error-2d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = 'V2',
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
   
end

function test_solve3d(nx, ny, nz, p, writeMatrix)
   local writeMatrix = writeMatrix or false
   print()
   print("Testing 3D Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {1.0, 1.0, 1.0},
      cells = {nx, ny, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 3, polyOrder = p }
   io.write("nx=",nx," ny=", ny," nz=", nz, " polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     bcLeft = { T = "D", V = 0.0 },
     bcRight = { T = "D", V = 0.0 },
     bcBottom = { T = "N", V = 0.0 },
     bcTop = { T = "D", V = 0.0 },
     -- periodicDirs = {1,2},
     writeStiffnessMatrix = writeMatrix,
   }
   local t2 = os.clock()
   io.write("3D Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1},
   }
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local y = xn[2]
                    local z = xn[3]
                    local a, b = 2, 5
                    local c1, d0 = 0, 0
                    local c0 = a/12 - 1/2
                    local d1 = b/12 - 1/2
                    local t1 = (1-a*x^2)*(-b*y^4/12 + y^2/2 + d0*y + d1)
                    local t2 = (1-b*y^2)*(-a*x^4/12 + x^2/2 + c0*x + c1)
                    return (z+1)*(t1+t2)
                 end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})
   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1},
   }
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x = xn[1]
                    local y = xn[2]
                    local z = xn[3]
                    local a, b = 2, 5
                    local c1, d0 = 0, 0
                    local c0 = a/12 - 1/2
                    local d1 = b/12 - 1/2
                    local t1 = x^2/2-a*x^4/12+c0*x+c1
                    local t2 = y^2/2-b*y^4/12+d0*y+d1
                    return (z+1)*t1*t2
                 end
   }
   initExactSolModal:advance(0.,0.,{},{exactSolModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1}
   }

   print("Solving...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("3D Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("phi-solution-3d.bp", 0.0)
   --exactSolModal:write("exact-solution-3d.bp", 0.0)
   --err:write("error-3d.bp", 0.0)
  
   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = 'V2',
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")

   return math.sqrt(lv[1])
end

function test_solve3d_periodic(nx, ny, nz, p)
   print()
   print("Testing 3D periodic Poisson solve...")
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {2*math.pi, 2*math.pi, 2*math.pi},
      cells = {nx, ny, nz},
   }
   local basis = Basis.CartModalSerendipity { ndim = 3, polyOrder = p }
   io.write("nx=",nx," ny=", ny, " nz=", nz, " polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemPerpPoisson {
     onGrid = grid,
     basis = basis,
     periodicDirs = {1,2},
     writeStiffnessMatrix = false,
   }
   local t2 = os.clock()
   io.write("3D periodic Poisson init took ", t2-t1, " s\n")
   local srcModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1},
   }
   -- initialization constants
   local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
   local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                 local x = xn[1]
                 local y = xn[2]
                 local t1, t2 = 0.0, 0.0
                 local f = 0.0
                 for m = 0,2 do
                    for n = 0,2 do
                       t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
                       t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
                       f = f + -(m*m+n*n)*(t1+t2)
                    end
                 end
                 return f/50.0
              end
   }
   initSrcModal:advance(0.,0.,{},{srcModal})

   -- calculate exact solution
   local exactSolModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1},
   }
   local initfunc = function (t,xn)
              local x = xn[1]
              local y = xn[2]
              local t1, t2 = 0.0, 0.0
              local f = 0.0
              for m = 0,2 do
                 for n = 0,2 do
                    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
                    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
                    f = f + (t1+t2)
                 end
              end
              f = f/50.0
              return f-.6
           end
   local initExactSolModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = initfunc,
   }
   initExactSolModal:advance(0.,0.,{},{exactSolModal})

   local phiModal = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1}
   }

   print("Solving...")
   local t1 = os.clock()
   poisson:advance(0.,0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("3D periodic Poisson solve took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1, 1}
   }

   err:combine(1.0, exactSolModal, -1.0, phiModal)

   --phiModal:write("periodic-phi-solution-3d.bp", 0.0)
   --exactSolModal:write("periodic-exact-solution-3d.bp", 0.0)
   --err:write("periodic-error-3d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
     onGrid = grid,
     basis = basis,
     numComponents = 1,
     quantity = 'V2',
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
   
end


function test_solve2d_p1()
  print("--- Testing convergence of 2D solver with p=1 ---")
  err1 = test_solve2d(32, 32, 1)
  err2 = test_solve2d(64, 64, 1)
  err3 = test_solve2d(128, 128, 1)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(1.0, err1/err2/4.0, .01)
  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

function test_solve2d_p2()
  print("--- Testing convergence of 2D solver with p=2 ---")
  err1 = test_solve2d(32, 32, 2)
  err2 = test_solve2d(64, 64, 2)
  err3 = test_solve2d(128, 128, 2)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(2.0, err1/err2/4.0, .01)
  assert_close(2.0, err2/err3/4.0, .01)
  print()
end

function test_solve3d_p1()
  print("--- Testing convergence of 3D solver with p=1 ---")
  err1 = test_solve3d(32, 32, 4, 1)
  err2 = test_solve3d(64, 64, 4, 1)
  err3 = test_solve3d(128, 128, 4, 1)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(1.0, err1/err2/4.0, .01)
  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

function test_solve3d_p2()
  print("--- Testing convergence of 3D solver with p=2 ---")
  err1 = test_solve3d(32, 32, 4, 2)
  err2 = test_solve3d(64, 64, 4, 2)
  err3 = test_solve3d(128, 128, 4, 2)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(2.0, err1/err2/4.0, .01)
  assert_close(2.0, err2/err3/4.0, .01)
  print()
end

function test_periodic2d_p1()
  print("--- Testing convergence of 2D periodic solver with p=1 ---")
  err1 = test_solve2d_periodic(32, 32, 1)
  err2 = test_solve2d_periodic(64, 64, 1)
  err3 = test_solve2d_periodic(128, 128, 1)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(1.0, err1/err2/4.0, .01)
  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

function test_periodic2d_p2()
  print("--- Testing convergence of 2D periodic solver with p=2 ---")
  err1 = test_solve2d_periodic(32, 32, 2)
  err2 = test_solve2d_periodic(64, 64, 2)
  err3 = test_solve2d_periodic(128, 128, 2)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(2.0, err1/err2/4.0, .01)
  assert_close(2.0, err2/err3/4.0, .01)
  print()
end

function test_periodic3d_p1()
  print("--- Testing convergence of 3D periodic solver with p=1 ---")
  err1 = test_solve3d_periodic(32, 32, 1, 1)
  err2 = test_solve3d_periodic(64, 64, 1, 1)
  err3 = test_solve3d_periodic(128, 128, 1, 1)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(1.0, err1/err2/4.0, .01)
  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

function test_periodic3d_p2()
  print("--- Testing convergence of 3D periodic solver with p=2 ---")
  err1 = test_solve3d_periodic(32, 32, 1, 2)
  err2 = test_solve3d_periodic(64, 64, 1, 2)
  err3 = test_solve3d_periodic(128, 128, 1, 2)
  print("Order:", err1/err2/4.0, err2/err3/4.0)
  assert_close(2.0, err1/err2/4.0, .01)
  assert_close(2.0, err2/err3/4.0, .01)
  print()
end

-- run tests
local t1 = os.clock()
test_init_nonperiodic()
test_init_allperiodic()
test_init_someperiodic()
test_init_errorcheck()
test_solve2d_p1()
test_solve2d_p2()
test_solve3d_p1()
test_solve3d_p2()
test_periodic2d_p1()
test_periodic2d_p2()
test_periodic3d_p1()
test_periodic3d_p2()
local t2 = os.clock()

print()
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
io.write("Total test time: ", t2-t1, " s\n")
