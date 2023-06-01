-- Gkyl ------------------------------------------------------------------------
--
-- Tests for FEM parallel Poisson solver
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Lin        = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
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

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData = {polyOrder  = basis:polyOrder(),
                  basisType  = basis:id(),},
   }
   return fld
end

function test_femparproj_3x(nx, ny, nz, p)
   print()
   print("Testing 3D parallel FEM projection...")

--   local lower = { 0.0,  0.0,  0.0}
--   local upper = { 1.0,  1.0,  1.0}
   local lower = {-2.0, -2.0, -.50}
   local upper = { 2.0,  2.0,  .50}
   local cells = {nx, ny, nz}
   local polyOrder = p

   local grid = createGrid(lower, upper, cells)
   local basis = createBasis(grid:ndim(), polyOrder)
   io.write("nx=",nx," ny=",ny," nz=",nz," polyOrder=", p, "\n")
   local t1 = os.clock()
   local poisson = Updater.FemParproj {
      onGrid  = grid, 
      basis   = basis,  
      periodicParallelDir = false,
   }
   local t2 = os.clock()
   local srcModal = createField(grid, basis)
   -- initialization constants
   local initSrcModal = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t,xn)
                    local x, y, z = xn[1], xn[2], xn[3]
                    --local a, b = 2, -12
                    --local c1 = 0
                    --local c0 = -a/12 - 1/2 - b/30
                    --return  1+a*z^2+b*z^4
		    local mu = {0.2, 0.2}
		    local sig = 0.3;
		    return math.exp(-(math.pow(x-mu[1],2)+math.pow(y-mu[2],2))/(2.0*sig*sig))*math.sin(2.*math.pi*z);
                 end
   }
   initSrcModal:advance(0.,{},{srcModal})

   local phiModal = createField(grid, basis)

   phiModal:copy(srcModal)
   print("Smoothing...")
   local t1 = os.clock()
   poisson:advance(0.,{srcModal},{phiModal})
   local t2 = os.clock()
   io.write("3D parallel Poisson smooth took total of ", t2-t1, " s\n")

   local err = createField(grid, basis)
   err:combine(1.0, srcModal, -1.0, phiModal)

--   srcModal:write("phi-rough-3d.bp", 0.0)
--   phiModal:write("phi-smooth-3d.bp", 0.0)
   --exactSolModal:write("exact-solution-1d.bp", 0.0)
   --err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   io.write("Average RMS error = ", math.sqrt(lv[1]), "\n")
   return math.sqrt(lv[1])
end

function test_3x_p1()
   print("--- Testing convergence of 3D FEM parallel projection with p=1 ---")
--   err1 = test_femparproj_3x(3, 3, 4, 1)
   err1 = test_femparproj_3x(4, 4, 32, 1)
   err2 = test_femparproj_3x(4, 4, 64, 1)
   err3 = test_femparproj_3x(4, 4, 128, 1)
   print("Order:", err1/err2/2.0, err2/err3/2.0)
   assert_close(4.0, err1/err2/2.0, .05)
   assert_close(4.0, err2/err3/2.0, .05)
   print()
end

local t1 = os.clock()
test_3x_p1()
local t2 = os.clock()
print()
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
io.write("Total test time: ", t2-t1, " s\n")
