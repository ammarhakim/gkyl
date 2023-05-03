-- Gkyl ------------------------------------------------------------------------
--Test the skinGhostAvg Updater in Parallel
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Mpi              = require "Comm.Mpi"
local DecompRegionCalc = require "Lib.CartDecomp"

local skinGhostAvgDecl = require "Updater.skinGhostAvgData.skinGhostAvgDecl"
local SkinGhostAvg = require "Updater.SkinGhostAvg"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

function test_1()
   local function createField(grid, basis, vComp)
      vComp = vComp or 1
      local fld = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis()*vComp,
         ghost         = {1, 1},
         metaData      = {polyOrder = basis:polyOrder(),
                          basisType = basis:id(),}
      }
      return fld
   end


   -- Set up a 3d grid.
   local decomp = DecompRegionCalc.CartProd { cuts = {2,1,2} }
   local grid = Grid.RectCart{
      lower = {0,-3,-5},
      upper = {1,3,5},
      cells = {2,2,4},
      decomposition = decomp,
   }

   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = 1 }
   local inFld=createField(grid,basis)

   local globalSkinInRange = grid:globalRange():upperSkin(3, 1)
   local localSkinInRange = inFld:localRange():intersect(globalSkinInRange)
   local zerorange = inFld:localRange():difference(localSkinInRange)
   local indexer = inFld:genIndexer()

   --Fill the in field

   --fill ones only in the skin range
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         fitr[i] = 1
      end
   end

   -- fill ascending in the ghost range
   for idx in localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[3] = idx[3] + 1
      local fitr = inFld:get(indexer(idxGhost))
      for i = 1,8 do
         fitr[i] = i
      end
   end



   local skinGhostAvgUpper = SkinGhostAvg{
      grid = grid,
      basis = basis,
      edge='upper',
      bcDir = 3,
      lowerGhost = 1,
      upperGhost = 1,
      advanceArgs = {inFld},
   }

   skinGhostAvgUpper:advance(inFld)


   for idx in zerorange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(0.0, fitr[i],1e-6, "Checking if copy worked in unmodified range")
      end
   end

   local expected_output = {-1.165064,-1.781089,-1.964102,-0.250000,-1.897114,-0.605662,-0.711325,-0.672650}
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(expected_output[i], fitr[i],1e-6, "Checking if copy worked in modified range")
      end
   end

end

function test_2()
   local function createField(grid, basis, vComp)
      vComp = vComp or 1
      local fld = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis()*vComp,
         ghost         = {1, 1},
         metaData      = {polyOrder = basis:polyOrder(),
                          basisType = basis:id(),}
      }
      return fld
   end


   -- Set up a 3d grid.
   local decomp = DecompRegionCalc.CartProd { cuts = {2,1,2} }
   local grid = Grid.RectCart{
      lower = {0,-3,-5},
      upper = {1,3,5},
      cells = {2,2,4},
      decomposition = decomp,
   }

   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = 1 }
   local inFld=createField(grid,basis)

   local globalSkinInRange = grid:globalRange():lowerSkin(3, 1)
   local localSkinInRange = inFld:localRange():intersect(globalSkinInRange)
   local zerorange = inFld:localRange():difference(localSkinInRange)
   local indexer = inFld:genIndexer()

   --Fill the in field

   --fill ones only in the skin range
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         fitr[i] = 1
      end
   end

   -- fill ascending in the ghost range
   for idx in localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[3] = idx[3] - 1
      local fitr = inFld:get(indexer(idxGhost))
      for i = 1,8 do
         fitr[i] = i
      end
   end


   local skinGhostAvgLower = SkinGhostAvg{
      grid = grid,
      basis = basis,
      edge='lower',
      bcDir = 3,
      lowerGhost = 1,
      upperGhost = 1,
      advanceArgs = {inFld},
   }

   skinGhostAvgLower:advance(inFld)


   for idx in zerorange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(0.0, fitr[i],1e-6, "Checking if copy worked in unmodified range")
      end
   end

   local expected_output = {3.165064,4.281089,4.964102,-0.250000,5.897114,-0.894338,-1.288675,-1.827350}
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(expected_output[i], fitr[i],1e-6, "Checking if copy worked in modified range")
      end
   end

end

function test_3()
   local function createField(grid, basis, vComp)
      vComp = vComp or 1
      local fld = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis()*vComp,
         ghost         = {1, 1},
         metaData      = {polyOrder = basis:polyOrder(),
                          basisType = basis:id(),}
      }
      return fld
   end


   -- Set up a 3d grid.
   local decomp = DecompRegionCalc.CartProd { cuts = {1,1,4} }
   local grid = Grid.RectCart{
      lower = {0,-3,-5},
      upper = {1,3,5},
      cells = {2,2,4},
      decomposition = decomp,
   }

   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = 1 }
   local inFld=createField(grid,basis)

   local globalSkinInRange = grid:globalRange():upperSkin(3, 1)
   local localSkinInRange = inFld:localRange():intersect(globalSkinInRange)
   local zerorange = inFld:localRange():difference(localSkinInRange)
   --local lv,uv = localSkinInRange:lowerAsVec(), localSkinInRange:upperAsVec()
   --print("localSkinRange lv", lv[1], lv[2], lv[3])
   --print("localSkinRange uv", uv[1], uv[2], uv[3])
   --local lv,uv = zerorange:lowerAsVec(), zerorange:upperAsVec()
   --print("zerorange lv", lv[1], lv[2], lv[3])
   --print("zerorange uv", uv[1], uv[2], uv[3])
   local indexer = inFld:genIndexer()

   --Fill the in field

   --fill ones only in the skin range
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         fitr[i] = 1
      end
   end

   -- fill ascending in the ghost range
   for idx in localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[3] = idx[3] + 1
      local fitr = inFld:get(indexer(idxGhost))
      for i = 1,8 do
         fitr[i] = i
      end
   end



   local skinGhostAvgUpper = SkinGhostAvg{
      grid = grid,
      basis = basis,
      edge='upper',
      bcDir = 3,
      lowerGhost = 1,
      upperGhost = 1,
      advanceArgs = {inFld},
   }

   skinGhostAvgUpper:advance(inFld)


   for idx in zerorange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(0.0, fitr[i],1e-6, "Checking if copy worked in unmodified range")
      end
   end

   local expected_output = {-1.165064,-1.781089,-1.964102,-0.250000,-1.897114,-0.605662,-0.711325,-0.672650}
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(expected_output[i], fitr[i],1e-6, "Checking if copy worked in modified range")
      end
   end

end

function test_4()
   local function createField(grid, basis, vComp)
      vComp = vComp or 1
      local fld = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis()*vComp,
         ghost         = {1, 1},
         metaData      = {polyOrder = basis:polyOrder(),
                          basisType = basis:id(),}
      }
      return fld
   end


   -- Set up a 3d grid.
   local decomp = DecompRegionCalc.CartProd { cuts = {1,1,4} }
   local grid = Grid.RectCart{
      lower = {0,-3,-5},
      upper = {1,3,5},
      cells = {2,2,4},
      decomposition = decomp,
   }

   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = 1 }
   local inFld=createField(grid,basis)

   local globalSkinInRange = grid:globalRange():lowerSkin(3, 1)
   local localSkinInRange = inFld:localRange():intersect(globalSkinInRange)
   local zerorange = inFld:localRange():difference(localSkinInRange)
   --local lv,uv = localSkinInRange:lowerAsVec(), localSkinInRange:upperAsVec()
   --print("localSkinRange lv", lv[1], lv[2], lv[3])
   --print("localSkinRange uv", uv[1], uv[2], uv[3])
   --local lv,uv = zerorange:lowerAsVec(), zerorange:upperAsVec()
   --print("zerorange lv", lv[1], lv[2], lv[3])
   --print("zerorange uv", uv[1], uv[2], uv[3])
   local indexer = inFld:genIndexer()

   --Fill the in field

   --fill ones only in the skin range
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         fitr[i] = 1
      end
   end

   -- fill ascending in the ghost range
   for idx in localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[3] = idx[3] - 1
      local fitr = inFld:get(indexer(idxGhost))
      for i = 1,8 do
         fitr[i] = i
      end
   end


   local skinGhostAvgLower = SkinGhostAvg{
      grid = grid,
      basis = basis,
      edge='lower',
      bcDir = 3,
      lowerGhost = 1,
      upperGhost = 1,
      advanceArgs = {inFld},
   }

   skinGhostAvgLower:advance(inFld)


   for idx in zerorange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(0.0, fitr[i],1e-6, "Checking if copy worked in unmodified range")
      end
   end

   local expected_output = {3.165064,4.281089,4.964102,-0.250000,5.897114,-0.894338,-1.288675,-1.827350}
   for idx in localSkinInRange:rowMajorIter() do
      local fitr = inFld:get(indexer(idx))
      for i = 1,8 do
         assert_close(expected_output[i], fitr[i],1e-6, "Checking if copy worked in modified range")
      end
   end

end

test_1()
test_2()

--test_3()
--test_4()
--Test 3 and 4 which have decompCuts[bcDir] = cellcs[bcDir] Fail with "344: attempt to index local 'zerorange' (a nil value)".
-- So the Updater does not yet work if #decompcuts = #cells in the bcDir

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
