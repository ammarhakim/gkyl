-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grid objects
--
-- See bottom of the file for a list of tests/select tests to run.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local Lin        = require "Lib.Linalg"
local Updater    = require "Updater"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local DecompRegionCalc = require "Lib.CartDecomp"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      cells = {10, 20}
   }

   -- Just make sure setIndex() method works. For RectCart object
   -- setting index is not needed.
   idx = Lin.IntVec(grid:ndim())
   idx[1], idx[2] = 1, 1
   grid:setIndex(idx)
   
   assert_equal(2, grid:ndim(), "Checking NDIM")
   assert_equal("uniform", grid:id(), "Checking ID")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")

   local localRange = grid:localRange()
   assert_equal(1, localRange:lower(1), "Checking region bounds")
   assert_equal(10, localRange:upper(1), "Checking region bounds")

   assert_equal(1, localRange:lower(2), "Checking region bounds")
   assert_equal(20, localRange:upper(2), "Checking region bounds")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx")
   assert_equal(0.05, grid:dx(2), "Checking dx")

   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)

   
   local xc = Lin.Vec(2)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)

	 assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	 assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
      end
   end
   
   grid:setPeriodicDirs({1})
   local periodicDirs = grid:getPeriodicDirs()
   assert_equal(periodicDirs[1], 1, "Checking periodic dirs")
   assert_equal(grid:isDirPeriodic(1), true, "Checking periodic dirs")
   assert_equal(grid:isDirPeriodic(2), false, "Checking periodic dirs")

   grid:setPeriodicDirs({2})
   local periodicDirs = grid:getPeriodicDirs()
   assert_equal(periodicDirs[1], 2, "Checking periodic dirs")
   assert_equal(grid:isDirPeriodic(1), false, "Checking periodic dirs")
   assert_equal(grid:isDirPeriodic(2), true, "Checking periodic dirs")
end

function test_2()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")

   local localRange = grid:localRange()
   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)

   local idx = Lin.IntVec(grid:ndim())
   local xc = Lin.Vec(2)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)

	 assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	 assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
      end
   end
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40}
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")

   local localRange = grid:localRange()
   assert_equal(1, localRange:lower(1), "Checking region bounds")
   assert_equal(10, localRange:upper(1), "Checking region bounds")

   assert_equal(1, localRange:lower(2), "Checking region bounds")
   assert_equal(20, localRange:upper(2), "Checking region bounds")

   assert_equal(1, localRange:lower(3), "Checking region bounds")
   assert_equal(40, localRange:upper(3), "Checking region bounds")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx = Lin.IntVec(grid:ndim())
   local xc = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
	    grid:setIndex( idx:setValues {i, j, k} )
	    grid:cellCenter(xc)
	    
	    assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	    assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
	    assert_equal(xc[3], loz+(k-0.5)*dz, "Testing cell-center coordinate")
	 end
      end
   end
end

function test_4()
   local grid = Grid.NonUniformRectCart {
      cells = {10, 10}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(10, grid:numCells(2), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx 1")
   assert_equal(0.1, grid:dx(2), "Checking dx 2")

   assert_equal(0.1*0.1, grid:cellVolume(), "Checking volume")

end

function test_5a()
   local grid = Grid.NonUniformRectCart {   
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40},
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates (THIS WORKS AS MESH IS ACTUALLY UNIFORM).
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx        = Lin.IntVec(grid:ndim())
   local localRange = grid:localRange()
   local xc         = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
            grid:setIndex( idx:setValues {i, j, k} )
            grid:cellCenter(xc)
        
            assert_equal(lox+(i-0.5)*dx, xc[1], "Testing cell-center x coordinate")
            assert_equal(loy+(j-0.5)*dy, xc[2], "Testing cell-center y coordinate")
            assert_equal(loz+(k-0.5)*dz, xc[3], "Testing cell-center z coordinate")
        
            local xcDir = {grid:cellCenterInDir(1), grid:cellCenterInDir(2), grid:cellCenterInDir(3)} 
            assert_equal(xc[1], xcDir[1], "Testing cell center in x direction")
            assert_equal(xc[2], xcDir[2], "Testing cell center in y direction")
            assert_equal(xc[3], xcDir[3], "Testing cell center in z direction")
            
            local loDir = {grid:cellLowerInDir(1), grid:cellLowerInDir(2), grid:cellLowerInDir(3)} 
            assert_equal(lox+(i-1)*dx, loDir[1], "Testing cell lower in x direction")
            assert_equal(loy+(j-1)*dy, loDir[2], "Testing cell lower in y direction")
            assert_equal(loz+(k-1)*dz, loDir[3], "Testing cell lower in z direction")
            
            local upDir = {grid:cellUpperInDir(1), grid:cellUpperInDir(2), grid:cellUpperInDir(3)} 
            assert_equal(lox+i*dx, upDir[1], "Testing cell upper in x direction")
            assert_equal(loy+j*dy, upDir[2], "Testing cell upper in y direction")
            assert_equal(loz+k*dz, upDir[3], "Testing cell upper in z direction")
	 end
      end
   end
end

function test_5()
   local grid = Grid.NonUniformRectCart {   
      lower    = {0.0, 1.0, 2.0},
      upper    = {2.0, 5.0, 10.0},
      cells    = {10, 20, 40},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return zeta
	 end,
	 function (zeta)
	    return zeta	    
	 end,
	 function (zeta)
	    return zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates (THIS WORKS AS MESH IS ACTUALLY UNIFORM).
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx        = Lin.IntVec(grid:ndim())
   local localRange = grid:localRange()
   local xc         = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
	    grid:setIndex( idx:setValues {i, j, k} )
	    grid:cellCenter(xc)

	    assert_equal(lox+(i-0.5)*dx, xc[1], "Testing cell-center x coordinate")
	    assert_equal(loy+(j-0.5)*dy, xc[2], "Testing cell-center y coordinate")
	    assert_equal(loz+(k-0.5)*dz, xc[3], "Testing cell-center z coordinate")

            local xcDir = {grid:cellCenterInDir(1), grid:cellCenterInDir(2), grid:cellCenterInDir(3)} 
	    assert_equal(xc[1], xcDir[1], "Testing cell center in x direction")
	    assert_equal(xc[2], xcDir[2], "Testing cell center in y direction")
	    assert_equal(xc[3], xcDir[3], "Testing cell center in z direction")
            
            local loDir = {grid:cellLowerInDir(1), grid:cellLowerInDir(2), grid:cellLowerInDir(3)} 
            assert_equal(lox+(i-1)*dx, loDir[1], "Testing cell lower in x direction")
            assert_equal(loy+(j-1)*dy, loDir[2], "Testing cell lower in y direction")
            assert_equal(loz+(k-1)*dz, loDir[3], "Testing cell lower in z direction")
            
            local upDir = {grid:cellUpperInDir(1), grid:cellUpperInDir(2), grid:cellUpperInDir(3)} 
            assert_equal(lox+i*dx, upDir[1], "Testing cell upper in x direction")
            assert_equal(loy+j*dy, upDir[2], "Testing cell upper in y direction")
            assert_equal(loz+k*dz, upDir[3], "Testing cell upper in z direction")
	 end
      end
   end   
end

function test_6()
   local grid = Grid.NonUniformRectCart {   
      cells    = {10, 20, 40},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return 2*zeta
	 end,
	 function (zeta)
	    return 2*zeta	    
	 end,
	 function (zeta)
	    return 2*zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")
   assert_equal(0.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(2.0, grid:upper(2), "Checking upper")
   assert_equal(2.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.1, grid:dx(2), "Checking dx")
   assert_equal(0.05, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.1*0.05, grid:cellVolume(), "Checking volume")
end

function test_7()
   local grid = Grid.NonUniformRectCart {   
      cells    = {3},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return zeta*zeta
	 end,
      }
   }

   local idx = Lin.IntVec(grid:ndim())
   
   assert_equal(1, grid:ndim(), "Checking NDIM")

   assert_equal(3, grid:numCells(1), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   grid:setIndex( idx:setValues {1} )
   assert_equal(1/9, grid:dx(1), "Checking dx")
   assert_equal(1/9, grid:cellVolume(), "Checking volume")

   grid:setIndex( idx:setValues {2} )
   assert_equal(1/3, grid:dx(1), "Checking dx")
   assert_equal(1/3, grid:cellVolume(), "Checking volume")

   grid:setIndex(idx:setValues {3} )
   assert_equal(5/9, grid:dx(1), "Checking dx")
   assert_equal(5/9, grid:cellVolume(), "Checking volume")
end

function test_8()
   local grid = Grid.RectCart {
      lower        = {0.0, 1.0, 1.0},
      upper        = {2.0, 5.0, 10.0},
      cells        = {10, 20, 30},
      periodicDirs = {1, 3},
   }

   -- Check periodicity.
   assert_equal(true, grid:isDirPeriodic(1), "Checking periodicity")
   assert_equal(false, grid:isDirPeriodic(2), "Checking periodicity")
   assert_equal(true, grid:isDirPeriodic(3), "Checking periodicity")
end

function test_9()
   local grid = Grid.NonUniformRectCart {
      lower        = {0.0, 1.0, 1.0},
      upper        = {2.0, 5.0, 10.0},
      cells        = {10, 20, 30},
      periodicDirs = {2, 3},
   }

   -- Check periodicity.
   assert_equal(false, grid:isDirPeriodic(1), "Checking periodicity (NU)")
   assert_equal(true, grid:isDirPeriodic(2), "Checking periodicity (NU)")
   assert_equal(true, grid:isDirPeriodic(3), "Checking periodicity (NU)")
end

function test_10()
   -- Test the creation of child grids from a 2D parent grid.
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1} }
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

   -- Test the creation of child grids from a 3D parent grid.
   local parentLower = {0.0, 1.0, -0.5}
   local parentUpper = {2.0, 5.0,  0.5}
   local parentCells = {10, 20, 16}
   local grid = Grid.RectCart {
      lower = parentLower,
      upper = parentUpper,
      cells = parentCells,
   }

   local childGrid = {}
   for d = 1, 3 do   -- Create an X sub-grid, a Y sub-grid and a Z sub-grid.
      local childGridIngr = grid:childGrid({d})
      childGrid[d] = Grid.RectCart {
         lower = childGridIngr.lower, 
         upper = childGridIngr.upper,
         cells = childGridIngr.cells,
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

   -- Create an YZ, XZ and XY subgrids.
   local childGrid2D = {}
   for cD = 1, 3 do
      local pD = {}
      for d = 1, 3 do pD[d] = d end
      table.remove(pD,cD)
      local childGridIngr = grid:childGrid(pD)
      childGrid2D[cD] = Grid.RectCart {
         lower = childGridIngr.lower, 
         upper = childGridIngr.upper,
         cells = childGridIngr.cells,
      }
      assert_equal(2, childGrid2D[cD]:ndim(), string.format("Checking 2D child %d NDIM",cD))
      for d = 1, 2 do
         assert_equal(parentCells[pD[d]], childGrid2D[cD]:numCells(d), string.format("Checking 2D child %d numCells",cD))
         assert_equal(parentLower[pD[d]], childGrid2D[cD]:lower(d), string.format("Checking 2D child %d lower",cD))
         assert_equal(parentUpper[pD[d]], childGrid2D[cD]:upper(d), string.format("Checking 2D child %d upper",cD))
         assert_equal((parentUpper[pD[d]]-parentLower[pD[d]])/parentCells[pD[d]], childGrid2D[cD]:dx(d), string.format("Checking 2D child %d dx",cD))
      end
      assert_equal( ((parentUpper[pD[1]]-parentLower[pD[1]])/parentCells[pD[1]])*((parentUpper[pD[2]]-parentLower[pD[2]])/parentCells[pD[2]]),
                    childGrid2D[cD]:cellVolume(), string.format("Checking 2D child %d volume",cD) )
      local childDim, localRange = childGrid2D[cD]:ndim(), childGrid2D[cD]:localRange()
      local lox, dx = childGrid2D[cD]:lower(1), childGrid2D[cD]:dx(1)
      local loy, dy = childGrid2D[cD]:lower(2), childGrid2D[cD]:dx(2)
      local idx, xc = Lin.IntVec(childDim), Lin.Vec(childDim)
      for i = localRange:lower(1), localRange:upper(1) do
         for j = localRange:lower(2), localRange:upper(2) do
            childGrid2D[cD]:setIndex( idx:setValues {i,j} )
            childGrid2D[cD]:cellCenter(xc)
            assert_equal(xc[1], lox+(i-0.5)*dx, string.format("Testing 2D child %d cell-center coordinate",cD))
	    assert_equal(xc[2], loy+(j-0.5)*dy, string.format("Testing 2D child %d cell-center coordinate",cD))
         end
      end
   end
end

function test_11()
   -- Test the grid's findCell method in 2D.
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   -- Find a point inside a cell.
   local fIdx, np = {0, 0}, {0.35, 2.1}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for interior point")
   grid:findCell(np, fIdx, true, {2,nil})    -- If we know fIdx[1].
   assert(2==fIdx[1] and 6==fIdx[2], "Checking 2D fIdx for interior point (known)")

   -- Find a point on x-boundary of two cells.
   local fIdx, np = {0, 0}, {0.4, 2.1}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for cell x-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for cell x-boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for cell x-boundary point (upper)")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for cell x-boundary point (upper)")
   grid:findCell(np, fIdx, nil, {nil,6})    -- If we know fIdx[2].
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for interior point (known)")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for interior point (known)")
   grid:findCell(np, fIdx, true, {2,nil})    -- If we know fIdx[1].
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for interior point (known)")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for interior point (known)")

   -- Find a point on y-boundary of two cells.
   local fIdx, np = {0, 0}, {0.35, 2.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for cell y-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 2D fIdx[2] for cell y-boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for cell y-boundary point (upper)")
   assert_equal(7, fIdx[2], "Checking 2D fIdx[2] for cell y-boundary point (upper)")

   -- Find a point on corner of 4 cells.
   local fIdx, np = {0, 0}, {0.4, 3.0}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for corner point (lower)")
   assert_equal(10, fIdx[2], "Checking 2D fIdx[2] for corner point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for corner point (upper)")
   assert_equal(11, fIdx[2], "Checking 2D fIdx[2] for corner point (upper)")

   -- Find a point on lower x domain boundary.
   local fIdx, np = {0, 0}, {0.0, 3.1}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 2D fIdx[1] for lower-x point")
   assert_equal(11, fIdx[2], "Checking 2D fIdx[2] for lower-x point")
   local fIdx, np = {0, 0}, {0.0, 3.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 2D fIdx[1] for lower-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking 2D fIdx[2] for lower-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(1, fIdx[1], "Checking 2D fIdx[1] for lower-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking 2D fIdx[2] for lower-x boundary point (upper)")

   -- Find a point on lower y domain boundary.
   local fIdx, np = {0, 0}, {0.5, 1.}
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for lower-y point")
   assert_equal(1, fIdx[2], "Checking 2D fIdx[2] for lower-y point")
   local fIdx, np = {0, 0}, {0.4, 1.}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for lower-y boundary point (lower)")
   assert_equal(1, fIdx[2], "Checking 2D fIdx[2] for lower-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for lower-y boundary point (upper)")
   assert_equal(1, fIdx[2], "Checking 2D fIdx[2] for lower-y boundary point (upper)")

   -- Find a point on upper x domain boundary.
   local fIdx, np = {0, 0}, {2.0, 3.1}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 2D fIdx[1] for upper-x point")
   assert_equal(11, fIdx[2], "Checking 2D fIdx[2] for upper-x point")
   local fIdx, np = {0, 0}, {2.0, 3.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 2D fIdx[1] for upper-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking 2D fIdx[2] for upper-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(10, fIdx[1], "Checking 2D fIdx[1] for upper-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking 2D fIdx[2] for upper-x boundary point (upper)")

   -- Find a point on upper y domain boundary.
   local fIdx, np = {0, 0}, {0.5, 5.}
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for upper-y point")
   assert_equal(20, fIdx[2], "Checking 2D fIdx[2] for upper-y point")
   local fIdx, np = {0, 0}, {0.4, 5.}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 2D fIdx[1] for upper-y boundary point (lower)")
   assert_equal(20, fIdx[2], "Checking 2D fIdx[2] for upper-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 2D fIdx[1] for upper-y boundary point (upper)")
   assert_equal(20, fIdx[2], "Checking 2D fIdx[2] for upper-y boundary point (upper)")

   -- Find domain corner points.
   local fIdx, np = {0, 0}, {0.0, 1.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower left corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower left corner point")
   local fIdx, np = {0, 0}, {2.0, 1.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for lower right corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower right corner point")
   local fIdx, np = {0, 0}, {0.0, 5.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for upper left corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper left corner point")
   local fIdx, np = {0, 0}, {2.0, 5.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper right corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper right corner point")

   -- Test the grid's findCell method in 3D.
   local grid = Grid.RectCart {
      lower = {0.0, 1.0, -0.5},
      upper = {2.0, 5.0,  0.5},
      cells = {10, 20, 8}
   }

   -- Find a point inside a cell.
   local fIdx, np = {0, 0, 0}, {0.35, 2.1, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for interior point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for interior point")
   grid:findCell(np, fIdx, nil, {2,nil,nil})   -- If we know fIdx[1].
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for interior point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for interior point")
   grid:findCell(np, fIdx, nil, {nil,6,nil})   -- If we know fIdx[2].
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for interior point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for interior point")
   grid:findCell(np, fIdx, nil, {nil,nil,3})   -- If we know fIdx[3].
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for interior point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for interior point")

   -- Find a point on x-boundary of two cells.
   local fIdx, np = {0, 0, 0}, {0.4, 2.1, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell x-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell x-boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell x-boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for cell x-boundary point (upper)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell x-boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell x-boundary point (upper)")
   grid:findCell(np, fIdx, nil, {2,nil,nil})   -- If we know fIdx[1].
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell x-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell x-boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell x-boundary point (lower)")
   grid:findCell(np, fIdx, false, {3,nil,nil})   -- If we know fIdx[1] and pickLower=False.
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for cell x-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell x-boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell x-boundary point (lower)")
   grid:findCell(np, fIdx, nil, {2,nil,3})   -- If we know fIdx[1] and fIdx[3].
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell x-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell x-boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell x-boundary point (lower)")

   -- Find a point on y-boundary of two cells.
   local fIdx, np = {0, 0, 0}, {0.35, 2.2, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell y-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell y-boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[2] for cell y-boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell y-boundary point (upper)")
   assert_equal(7, fIdx[2], "Checking 3D fIdx[2] for cell y-boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for cell y-boundary point (upper)")

   -- Find a point on z-boundary of two cells.
   local fIdx, np = {0, 0, 0}, {0.35, 2.1, 0.}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell z-boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell z-boundary point (lower)")
   assert_equal(4, fIdx[3], "Checking 3D fIdx[2] for cell z-boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for cell z-boundary point (upper)")
   assert_equal(6, fIdx[2], "Checking 3D fIdx[2] for cell z-boundary point (upper)")
   assert_equal(5, fIdx[3], "Checking 3D fIdx[3] for cell z-boundary point (upper)")

   -- Find a point on corner of 4 cells.
   local fIdx, np = {0, 0, 0}, {0.4, 3.0, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for 4-corner point (lower)")
   assert_equal(10, fIdx[2], "Checking 3D fIdx[2] for 4-corner point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for 4-corner point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for 4-corner point (upper)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for 4-corner point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for 4-corner point (upper)")

   -- Find a point on corner of 8 cells.
   local fIdx, np = {0, 0, 0}, {0.4, 3.0, 0.125}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for 8-corner point (lower)")
   assert_equal(10, fIdx[2], "Checking 3D fIdx[2] for 8-corner point (lower)")
   assert_equal(5, fIdx[3], "Checking 3D fIdx[3] for 8-corner point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for 8-corner point (upper)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for 8-corner point (upper)")
   assert_equal(6, fIdx[3], "Checking 3D fIdx[3] for 8-corner point (upper)")

   -- Find a point on lower x domain boundary.
   local fIdx, np = {0, 0, 0}, {0.0, 3.1, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower-x point")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-x point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-x point")
   local fIdx, np = {0, 0, 0}, {0.0, 3., -0.2}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking 3D fIdx[2] for lower-x boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-x boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-x boundary point (upper)")

   -- Find a point on lower y domain boundary.
   local fIdx, np = {0, 0, 0}, {0.5, 1., -0.2}
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-y point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower-y point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-y point")
   local fIdx, np = {0, 0, 0}, {0.4, 1., -0.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for lower-y boundary point (lower)")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower-y boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-y boundary point (upper)")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower-y boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for lower-y boundary point (upper)")

   -- Find a point on lower z domain boundary.
   local fIdx, np = {0, 0, 0}, {0.5, 3.1, -0.5} -- Inside a cell.
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-z point")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-z point")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower-z point")
   local fIdx, np = {0, 0, 0}, {0.4, 3.1, -0.5}  -- On x-boundary of 2 cells.
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for lower-z boundary point (lower-x)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-z boundary point (lower-x)")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower-z boundary point (lower-x)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-z boundary point (upper-x)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-z boundary point (upper-x)")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower-z boundary point (upper-x)")
   local fIdx, np = {0, 0, 0}, {0.5, 3., -0.5}  -- On y-boundary of 2 cells.
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-z boundary point (lower-y)")
   assert_equal(10, fIdx[2], "Checking 3D fIdx[2] for lower-z boundary point (lower-y)")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower-z boundary point (lower-y)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for lower-z boundary point (upper-y)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for lower-z boundary point (upper-y)")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower-z boundary point (upper-y)")

   -- Find a point on upper x domain boundary.
   local fIdx, np = {0, 0, 0}, {2.0, 3.1, -0.2}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper-x point")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for upper-x point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-x point")
   local fIdx, np = {0, 0, 0}, {2.0, 3., -0.2}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking 3D fIdx[2] for upper-x boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for upper-x boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-x boundary point (upper)")

   -- Find a point on upper y domain boundary.
   local fIdx, np = {0, 0, 0}, {0.5, 5., -0.2}   -- Inside a cell.
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-y point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper-y point")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-y point")
   local fIdx, np = {0, 0, 0}, {0.4, 5., -0.2}   -- On x-boundary of 2 cells.
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for upper-y boundary point (lower)")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper-y boundary point (lower)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-y boundary point (upper)")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper-y boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-y boundary point (upper)")
   local fIdx, np = {0, 0, 0}, {0.5, 5., -0.25}   -- On z-boundary of 2 cells.
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-y boundary point (lower)")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper-y boundary point (lower)")
   assert_equal(2, fIdx[3], "Checking 3D fIdx[3] for upper-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-y boundary point (upper)")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper-y boundary point (upper)")
   assert_equal(3, fIdx[3], "Checking 3D fIdx[3] for upper-y boundary point (upper)")

   -- Find a point on upper z domain boundary.
   local fIdx, np = {0, 0, 0}, {0.5, 3.1, 0.5}   -- Inside a cell.
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-z point")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for upper-z point")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for upper-z point")
   local fIdx, np = {0, 0, 0}, {0.4, 3.1, 0.5}   -- On x-boundary of 2 cells.
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking 3D fIdx[1] for upper-z boundary point (lower)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for upper-z boundary point (lower)")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for upper-z boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking 3D fIdx[1] for upper-z boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking 3D fIdx[2] for upper-z boundary point (upper)")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for upper-z boundary point (upper)")

   -- Find domain corner points.
   local fIdx, np = {0, 0, 0}, {0.0, 1., -0.5}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower left bottom corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower left bottom corner point")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower left bottom corner point")
   local fIdx, np = {0, 0, 0}, {2.0, 1., -0.5}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for lower right bottom corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower right bottom corner point")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for lower right bottom corner point")
   local fIdx, np = {0, 0, 0}, {0.0, 5., -0.5}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for upper left bottom corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper left bottom corner point")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for upper left bottom corner point")
   local fIdx, np = {0, 0, 0}, {2.0, 5., -0.5}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper right bottom corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper right bottom corner point")
   assert_equal(1, fIdx[3], "Checking 3D fIdx[3] for upper right bottom corner point")
   local fIdx, np = {0, 0, 0}, {0.0, 1., 0.5}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for lower left top corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower left top corner point")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for lower left top corner point")
   local fIdx, np = {0, 0, 0}, {2.0, 1., 0.5}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for lower right top corner point")
   assert_equal(1, fIdx[2], "Checking 3D fIdx[2] for lower right top corner point")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for lower right top corner point")
   local fIdx, np = {0, 0, 0}, {0.0, 5., 0.5}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking 3D fIdx[1] for upper left top corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper left top corner point")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for upper left top corner point")
   local fIdx, np = {0, 0, 0}, {2.0, 5., 0.5}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking 3D fIdx[1] for upper right top corner point")
   assert_equal(20, fIdx[2], "Checking 3D fIdx[2] for upper right top corner point")
   assert_equal(8, fIdx[3], "Checking 3D fIdx[3] for upper right top corner point")
end

-- Run tests.
test_1()
test_2()
test_3()
test_4()
test_5a()
test_5()
test_6()
test_7()
test_8()
test_9()
test_10()
test_11()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
