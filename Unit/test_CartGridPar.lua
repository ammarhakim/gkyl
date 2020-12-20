-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grids in parallel
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local Lin        = require "Lib.Linalg"
local Mpi        = require "Comm.Mpi"
local Updater    = require "Updater"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local DecompRegionCalc = require "Lib.CartDecomp"

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
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {1.0, 1.0},
      cells = {100, 200},
      decomposition = decomp,
   }

   local localRange = grid:localRange()
   assert_equal(50*100, localRange:volume(), "Checking volume of local ranges")
   assert_equal(rnk+1, grid:subGridId(), "Checking sub-grid ID is correct")
   assert_equal(false, grid:isShared(), "Checking if grid is shared")
end

function test_2(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1}, useShared = true }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {1.0, 1.0},
      cells = {100, 200},
      decomposition = decomp,
   }

   local localRange = grid:localRange()
   assert_equal(100*200, localRange:volume(), "Checking volume of local ranges")
   assert_equal(1, grid:subGridId(), "Checking sub-grid ID is correct")
   assert_equal(true, grid:isShared(), "Checking if grid is shared")   
end

function test_3(comm)
   -- Test the creation of child grids from a 2D parent grid.
   local sz, rnk  = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   local parentLower = {0.0, 1.0}
   local parentUpper = {2.0, 5.0}
   local parentCells = {10, 20}
   local grid = Grid.RectCart {
      lower = parentLower,
      upper = parentUpper,
      cells = parentCells,
      decomposition = decomp,
   }

   local childGrid = {}
   for d = 1, 2 do   -- Create an X sub-grid and a Y sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower,
         upper = childGridIngr.upper,
         cells = childGridIngr.cells,
         decomposition = childGridIngr.decomposition,
      }
      assert_equal(1, childGrid[d]:ndim(), string.format("Checking child %d NDIM",d))
      assert_equal(parentCells[d], childGrid[d]:numCells(1), string.format("Checking child %d numCells",d))
      assert_equal(parentLower[d], childGrid[d]:lower(1), string.format("Checking child %d lower",d))
      assert_equal(parentUpper[d], childGrid[d]:upper(1), string.format("Checking child %d upper",d))
      assert_equal((parentUpper[d]-parentLower[d])/parentCells[d], childGrid[d]:dx(1), string.format("Checking child %d dx",d))
      assert_equal((parentUpper[d]-parentLower[d])/parentCells[d], childGrid[d]:cellVolume(), string.format("Checking child %d volume",d))
      local childDim, localRange = childGrid[d]:ndim(), childGrid[d]:localRange()
      local lox, dx = childGrid[d]:lower(1), childGrid[d]:dx(1)
      local idx, xc = Lin.IntVec(childDim), Lin.Vec(childDim)
      for i = localRange:lower(1), localRange:upper(1) do
         childGrid[d]:setIndex( idx:setValues {i} )
         childGrid[d]:cellCenter(xc)
         assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing child %d cell-center coordinate",d))
      end
   end

   -- Project a function on the 2D grid, then project onto a 1D grid.
   local polyOrder = 2
   local testFuncDir = {}
   testFuncDir[1] = function(x) return math.cos((2.*math.pi/2.)*x) end
   testFuncDir[2] = function(y) return math.sin((2.*math.pi/6.)*y) end
   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
   local testFunc  = function(t, xn)
      local x, y = xn[1], xn[2]
      return testFuncDir[1](x)*testFuncDir[2](y)
   end
   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()}
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = testFunc
   }
   project:advance(0., {}, {distf})
--   distf:write("distf.bp")
   local childBasis, distfChild, testFuncChild, projectChild = {}, {}, {}, {}
   for d = 1, 2 do
      childBasis[d] = Basis.CartModalSerendipity { ndim = childGrid[d]:ndim(), polyOrder = polyOrder }
      distfChild[d] = DataStruct.Field {
         onGrid        = childGrid[d],
         numComponents = childBasis[d]:numBasis(),
         ghost         = {1, 1},
         metaData      = {polyOrder = childBasis[d]:polyOrder(), basisType = childBasis[d]:id()}
      }
      testFuncChild[d] = function(t, xn) return testFuncDir[d](xn[1]) end
      projectChild[d] = Updater.ProjectOnBasis {
         onGrid   = childGrid[d],
         basis    = childBasis[d],
         evaluate = testFuncChild[d]
      }
      projectChild[d]:advance(0., {}, {distfChild[d]})
--      distfChild[d]:write("distfChild-" .. d ..".bp")  -- Note: if writing different datasets, cannot use _ to differentiate files.
   end

   Mpi.Barrier(comm)
end

-- Run tests
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)

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
