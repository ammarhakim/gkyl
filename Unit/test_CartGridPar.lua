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
      lower = parentLower,  cells         = parentCells,
      upper = parentUpper,  decomposition = decomp,
   }

   local childGrid = {}
   for d = 1, 2 do   -- Create an X sub-grid and a Y sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
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
      -- MF 2021/09/15 WARNING: this test is not sufficient. It doesn't check that the range is correct.
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

function test_4(comm)
   -- Test the creation of child grids from a 3D parent grid.
   local sz, rnk  = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1, 2} }
   local parentLower = {0.0, 1.0, -0.5}
   local parentUpper = {2.0, 5.0,  0.5}
   local parentCells = {10, 20, 16}
   local grid = Grid.RectCart {
      lower = parentLower,  cells         = parentCells,
      upper = parentUpper,  decomposition = decomp,
   }
   local ndim = #parentCells

   local childGrid = {}
   for d = 1, ndim do   -- Create an X sub-grid, a Y sub-grid and a Z sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
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
      -- MF 2021/09/15 WARNING: this test is not sufficient. It doesn't check that the range is correct.
      for i = localRange:lower(1), localRange:upper(1) do
         childGrid[d]:setIndex( idx:setValues {i} )
         childGrid[d]:cellCenter(xc)
         assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing child %d cell-center coordinate",d))
      end
   end

   -- Create an YZ, XZ and XY subgrids.
   local childGrid2D = {}
   for cD = 1, ndim do
      local pD = {}
      for d = 1, ndim do pD[d] = d end
      table.remove(pD,cD)
      local childGridIngr = grid:childGrid(pD)
      childGrid2D[cD] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
      }
      assert_equal(ndim-1, childGrid2D[cD]:ndim(), string.format("Checking 2D child %d NDIM",cD))
      for d = 1, ndim-1 do
         assert_equal(parentCells[pD[d]], childGrid2D[cD]:numCells(d), string.format("Checking 2D child %d numCells",cD))
         assert_equal(parentLower[pD[d]], childGrid2D[cD]:lower(d), string.format("Checking 2D child %d lower",cD))
         assert_equal(parentUpper[pD[d]], childGrid2D[cD]:upper(d), string.format("Checking 2D child %d upper",cD))
         assert_equal((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]], childGrid2D[cD]:dx(d), string.format("Checking 2D child %d dx",cD))
      end
      local vol = 1
      for d = 1, ndim-1 do vol = vol*((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]]) end
      assert_equal( vol, childGrid2D[cD]:cellVolume(), string.format("Checking 2D child %d volume",cD) )
      local childDim, localRange = childGrid2D[cD]:ndim(), childGrid2D[cD]:localRange()
      local lox, dx = childGrid2D[cD]:lower(1), childGrid2D[cD]:dx(1)
      local loy, dy = childGrid2D[cD]:lower(2), childGrid2D[cD]:dx(2)
      local idx, xc = Lin.IntVec(childDim), Lin.Vec(childDim)
      -- MF 2021/09/15 WARNING: this test is not sufficient. It doesn't check that the range is correct.
      for i = localRange:lower(1), localRange:upper(1) do
         for j = localRange:lower(2), localRange:upper(2) do
            childGrid2D[cD]:setIndex( idx:setValues {i,j} )
            childGrid2D[cD]:cellCenter(xc)
            assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing 2D child %d cell-center coordinate",cD))
	    assert_equal(xc[2], loy+(j-0.5)*dy, string.format("Testing 2D child %d cell-center coordinate",cD))
         end
      end
   end

   Mpi.Barrier(comm)
end

function test_5(comm)
   -- Test the creation of child grids from a 4D parent grid.
   local sz, rnk  = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1, 2, 1} }
   local parentLower = {0.0, 1.0, -0.5, -1.0}
   local parentUpper = {2.0, 5.0,  0.5,  1.0}
   local parentCells = {10, 20, 16, 12}
   local grid = Grid.RectCart {
      lower = parentLower,  cells         = parentCells,
      upper = parentUpper,  decomposition = decomp,
   }
   local ndim = #parentCells

   local childGrid = {}
   for d = 1, ndim do   -- Create an X, a Y, a Z and a Vx sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
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

   -- Create an YZVx, XZVx, XYVx and XYZ subgrids.
   local childGrid3D = {}
   for cD = 1, ndim do
      local pD = {}
      for d = 1, ndim do pD[d] = d end
      table.remove(pD,cD)
      local childGridIngr = grid:childGrid(pD)
      childGrid3D[cD] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
      }
      assert_equal(ndim-1, childGrid3D[cD]:ndim(), string.format("Checking 3D child %d NDIM",cD))
      for d = 1, ndim-1 do
         assert_equal(parentCells[pD[d]], childGrid3D[cD]:numCells(d), string.format("Checking 3D child %d numCells",cD))
         assert_equal(parentLower[pD[d]], childGrid3D[cD]:lower(d), string.format("Checking 3D child %d lower",cD))
         assert_equal(parentUpper[pD[d]], childGrid3D[cD]:upper(d), string.format("Checking 3D child %d upper",cD))
         assert_equal((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]], childGrid3D[cD]:dx(d), string.format("Checking 3D child %d dx",cD))
      end
      local vol = 1
      for d = 1, ndim-1 do vol = vol*((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]]) end
      assert_equal( vol, childGrid3D[cD]:cellVolume(), string.format("Checking 3D child %d volume",cD) )
      local childDim, localRange = childGrid3D[cD]:ndim(), childGrid3D[cD]:localRange()
      local lox, dx = childGrid3D[cD]:lower(1), childGrid3D[cD]:dx(1)
      local loy, dy = childGrid3D[cD]:lower(2), childGrid3D[cD]:dx(2)
      local loz, dz = childGrid3D[cD]:lower(3), childGrid3D[cD]:dx(3)
      local idx, xc = Lin.IntVec(childDim), Lin.Vec(childDim)
      -- MF 2021/09/15 WARNING: this test is not sufficient. It doesn't check that the range is correct.
      for i = localRange:lower(1), localRange:upper(1) do
         for j = localRange:lower(2), localRange:upper(2) do
            for k = localRange:lower(3), localRange:upper(3) do
               childGrid3D[cD]:setIndex( idx:setValues {i,j,k} )
               childGrid3D[cD]:cellCenter(xc)
               assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing 3D child %d cell-center coordinate",cD))
	       assert_equal(xc[2], loy+(j-0.5)*dy, string.format("Testing 3D child %d cell-center coordinate",cD))
	       assert_equal(xc[3], loz+(k-0.5)*dz, string.format("Testing 3D child %d cell-center coordinate",cD))
            end
         end
      end
   end

   Mpi.Barrier(comm)
end

function test_6(comm)
   -- Test the creation of child grids from a 5D parent grid.
   local sz, rnk  = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for CartGrid not run as number of procs not exactly 4")
      return
   end
   
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1, 1, 1, 1} }
   local parentLower = {0.0, 1.0, -0.5, -1.0, 0.}
   local parentUpper = {2.0, 5.0,  0.5,  1.0, 3.}
   local parentCells = {10, 20, 16, 12, 6}
   local grid = Grid.RectCart {
      lower = parentLower,  cells         = parentCells,
      upper = parentUpper,  decomposition = decomp,
   }
   local ndim = #parentCells

   local childGrid = {}
   for d = 1, ndim do   -- Create an X, a Y, a Z, a Vx and a Vy sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
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

   -- Create an YZVxVy, XZVxVy, XYVxVy, XYZVy and XYZVx subgrids.
   local childGrid4D = {}
   for cD = 1, ndim do
      local pD = {}
      for d = 1, ndim do pD[d] = d end
      table.remove(pD,cD)
      local childGridIngr = grid:childGrid(pD)
      childGrid4D[cD] = Grid.RectCart {
         lower = childGridIngr.lower,  cells         = childGridIngr.cells,
         upper = childGridIngr.upper,  decomposition = childGridIngr.decomposition,
      }
      assert_equal(ndim-1, childGrid4D[cD]:ndim(), string.format("Checking 4D child %d NDIM",cD))
      for d = 1, ndim-1 do
         assert_equal(parentCells[pD[d]], childGrid4D[cD]:numCells(d), string.format("Checking 4D child %d numCells",cD))
         assert_equal(parentLower[pD[d]], childGrid4D[cD]:lower(d), string.format("Checking 4D child %d lower",cD))
         assert_equal(parentUpper[pD[d]], childGrid4D[cD]:upper(d), string.format("Checking 4D child %d upper",cD))
         assert_equal((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]], childGrid4D[cD]:dx(d), string.format("Checking 4D child %d dx",cD))
      end
      local vol = 1
      for d = 1, ndim-1 do vol = vol*((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]]) end
      assert_equal( vol, childGrid4D[cD]:cellVolume(), string.format("Checking 4D child %d volume",cD) )
      local childDim, localRange = childGrid4D[cD]:ndim(), childGrid4D[cD]:localRange()
      local lox, dx = childGrid4D[cD]:lower(1), childGrid4D[cD]:dx(1)
      local loy, dy = childGrid4D[cD]:lower(2), childGrid4D[cD]:dx(2)
      local loz, dz = childGrid4D[cD]:lower(3), childGrid4D[cD]:dx(3)
      local lovx, dvx = childGrid4D[cD]:lower(4), childGrid4D[cD]:dx(4)
      local idx, xc = Lin.IntVec(childDim), Lin.Vec(childDim)
      -- MF 2021/09/15 WARNING: this test is not sufficient. It doesn't check that the range is correct.
      for i = localRange:lower(1), localRange:upper(1) do
         for j = localRange:lower(2), localRange:upper(2) do
            for k = localRange:lower(3), localRange:upper(3) do
               for l = localRange:lower(3), localRange:upper(3) do
                  childGrid4D[cD]:setIndex( idx:setValues {i,j,k,l} )
                  childGrid4D[cD]:cellCenter(xc)
                  assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing 4D child %d cell-center coordinate",cD))
	          assert_equal(xc[2], loy+(j-0.5)*dy, string.format("Testing 4D child %d cell-center coordinate",cD))
	          assert_equal(xc[3], loz+(k-0.5)*dz, string.format("Testing 4D child %d cell-center coordinate",cD))
	          assert_equal(xc[4], lovx+(l-0.5)*dvx, string.format("Testing 4D child %d cell-center coordinate",cD))
               end
            end
         end
      end
   end

   Mpi.Barrier(comm)
end

function test_7(comm)
   -- Test the grid's findCell method in 2D.

   local nRanks = Mpi.Comm_size(Mpi.COMM_WORLD)
   if nRanks ~= 2 then
      log("Not running test_15 as numProcs not exactly 4")
      return
   end

   local decomp = DecompRegionCalc.CartProd { cuts = {2, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20},
      decomposition = decomp,
   }

   -- Find a point inside a cell.
   local fIdx, np = {0, 0}, {0.35, 2.1}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for interior point")
   grid:findCell(np, fIdx, true, {2,nil})    -- If we know fIdx[1].

   assert(2==fIdx[1] and 6==fIdx[2], "Checking 2D fIdx for interior point (known)")
end

-- Run tests
--test_1(Mpi.COMM_WORLD)
--test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)
test_4(Mpi.COMM_WORLD)
test_5(Mpi.COMM_WORLD)
test_6(Mpi.COMM_WORLD)
test_7(Mpi.COMM_WORLD)

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
