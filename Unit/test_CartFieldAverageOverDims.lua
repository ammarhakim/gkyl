-- Gkyl ------------------------------------------------------------------------
--
-- Test updater that computes averages of a field over some dimensions,
-- CartFieldAverageOverDims.
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
local Time       = require "Lib.Time"

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
   elseif (bKind=="Max") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
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
      metaData = {
         polyOrder  = basis:polyOrder(),
         basisType  = basis:id(),
      }
   }
   return fld
end

function test_2x()
   local lower    = {-math.pi, -math.pi}
   local upper    = { math.pi,  math.pi}
   local numCells = {6, 8}

   -- Tolerances for polyOrder=1-3.
   local tols = {{1.e-10,1.e-10},{1.e-7,1.e-7},{1.e-7,1.e-6}}

   for polyOrder = 1, 2 do

      local grid  = createGrid(lower, upper, numCells)
      local basis = createBasis(grid:ndim(), polyOrder)
      local fld = createField(grid, basis)
      local projScalarFunc = Updater.ProjectOnBasis {
         onGrid   = grid,
         basis    = basis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local fFunc = function(t, xn)  -- Function to project onto basis.
         local x, y = xn[1], xn[2]
         return (1.+0.5*math.sin(x))*(1.+math.cos(y))
      end
      projScalarFunc:setFunc(fFunc)
      projScalarFunc:advance(0., {}, {fld})

      for avgDir = 1, 2 do

         -- Create child grid, child field and project expected answer.
         local gridIngr = grid:childGrid({avgDir==1 and 2 or 1})
         local childGrid = createGrid(gridIngr.lower, gridIngr.upper, gridIngr.cells)
         local childBasis = createBasis(childGrid:ndim(), polyOrder)
         local avgFldProj = createField(childGrid, childBasis)
         local childProjScalarFunc = Updater.ProjectOnBasis {
            onGrid   = childGrid,
            basis    = childBasis,
            evaluate = function (t, xn) return 1.0 end   -- Set later.
         }
         local avgfFunc = avgDir == 1
            and function(t, xn)  -- Function to project onto basis.
               local y = xn[1]
               return 1.+math.cos(y)
            end
            or function(t, xn)  -- Function to project onto basis.
               local x = xn[1]
               return 1.+0.5*math.sin(x)
            end
         childProjScalarFunc:setFunc(avgfFunc)
         childProjScalarFunc:advance(0., {}, {avgFldProj})

         local avgFld = createField(childGrid, childBasis)
   
         local avgUpd = Updater.CartFieldAverageOverDims {
            onGrid      = grid,       onBasis    = basis,
            childGrid   = childGrid,  childBasis = childBasis,
            averageDirs = {avgDir},
         }
   
         avgUpd:advance(0., {fld}, {avgFld})  -- Compute the averaged f.
   
         local indexer = avgFld:genIndexer()
         local localRange = avgFld:localRange()
         local avgFldPtr, avgFldProjPtr = avgFld:get(1), avgFldProj:get(1)
         for idx in localRange:rowMajorIter() do
            avgFld:fill(indexer(idx), avgFldPtr)
            avgFldProj:fill(indexer(idx), avgFldProjPtr)
            for k = 1, childBasis:numBasis() do
               assert_close(avgFldProjPtr[1], avgFldPtr[1], tols[polyOrder][1], string.format("Checking 2x p%d avg(f) vs ProjectOnBasis.",polyOrder))
            end
         end
      end
   end
end

-- Run tests.
test_2x()
--test_3x()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
