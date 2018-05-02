-- Gkyl ------------------------------------------------------------------------
--
-- Tests for computing additional moments needed by Lenard-Bernstein operator. 
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

function test_moms1x1v(nx, nv, p, writeMatrix)
   writeMatrix = writeMatrix or false
   print()
   print("Testing 1x1v moment calculation...")
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {nx, nv},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nx=",nx," nv=",nv," polyOrder=", p, "\n")

   -- Zeroth and first (parallel) moments.
   local numDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1par = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local mom2gk     = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Calculated and analytica parallel flow speed.
   local Upar = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local UparA = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Initialize number density and first (parallel) moment.
   local initGkDens = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return 0.90-x
                    return math.exp(-x)-0.5
                 end
   }
   initGkDens:advance(0.,0.,{},{numDens})
   local initGkUpar = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                      return (math.exp(-((x-mu)/(math.sqrt(2)*sig))^2))*math.tanh(x-0.5)
--                      return (0.90-x)*math.tanh(x-0.5)
--                    return (0.90-x)*math.sin(4.0*2.0*math.pi*x)
                    return (math.exp(-x)-0.5)*math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initGkUpar:advance(0.,0.,{},{mom1par})
   -- Analytic parallel flow speed.
   local initUpar = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return math.tanh(x-0.5)
                    return math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initUpar:advance(0.,0.,{},{UparA})

   -- Compute parallel frow speed.
   local calcUpar = Updater.DistFuncPrimMomentCalc {
      onGrid     = confGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "GkUpar",
   }
--   initUpar:advance(0.,0.,{},{Upar})
   calcUpar:advance(0.,0.,{numDens,mom1par,mom2gk},{Upar})

   print("Computing parallel flow speed...")
   local t1 = os.clock()
   local t2 = os.clock()
   io.write("Upar computation took total of ", t2-t1, " s\n")

   local err = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost = {1, 1}
   }

   err:combine(1.0, UparA, -1.0, Upar)

   numDens:write("numDens-1d.bp", 0.0)
   mom1par:write("mom1par-1d.bp", 0.0)
   Upar:write("Upar-1d.bp", 0.0)
   UparA:write("UparA-1d.bp", 0.0)
   --err:write("error-2d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = confBasis,
      numComponents = 1,
      quantity      = "V2"
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_moms1x1v_p1()
  print("--- Testing convergence of 2D solver with p=1 ---")
  err1 = test_moms1x1v(32, 8, 1)
--  err2 = test_moms1x1v(64, 1)
--  err3 = test_moms1x1v(128, 1)
--  print("Order:", err1/err2/4.0, err2/err3/4.0)
--  assert_close(1.0, err1/err2/4.0, .01)
--  assert_close(1.0, err2/err3/4.0, .01)
--  print()
end

-- function test_solve2d_p2()
--   print("--- Testing convergence of 2D solver with p=2 ---")
--   err1 = test_solve2d(32, 32, 2)
--   err2 = test_solve2d(64, 64, 2)
--   err3 = test_solve2d(128, 128, 2)
--   print("Order:", err1/err2/4.0, err2/err3/4.0)
--   assert_close(2.0, err1/err2/4.0, .01)
--   assert_close(2.0, err2/err3/4.0, .01)
--   print()
-- end

-- run tests
local t1 = os.clock()
test_moms1x1v_p1()
-- test_solve2d_p2()
local t2 = os.clock()

print()
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
io.write("Total test time: ", t2-t1, " s\n")
