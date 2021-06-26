
local Unit       = require "Unit"
local Basis      = require "Basis"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test()
   local grid = Grid.RectCart {
      lower = {0.0, },
      upper = {1.0, },
      cells = {  2, },
   }
   local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = 2 }
   local numVals = 2
   local vecField = DataStruct.Field {
      onGrid        = grid,
      numComponents = numVals*basis:numBasis(),
      ghost         = {0, 0},
   }
   local compFields = {}
   for i = 1, numVals do
      compFields[i] = DataStruct.Field {
         onGrid        = grid,
         numComponents = basis:numBasis(),
         ghost         = {0, 0},
      }
   end
   
   local localRange = vecField:localRange()
   local indexer    = vecField:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = vecField:get(indexer(idx))
      for i = 1, numVals*basis:numBasis() do
         fitr[i] = idx[1]+i
      end
   end

   local separate = Updater.SeparateVectorComponents {
      onGrid = grid,
      basis  = basis,
   }

   separate:advance(0, {vecField}, compFields)

   local cIndexer = compFields[1]:genIndexer()
   for idx in localRange:colMajorIter() do
      local vItr = vecField:get(indexer(idx))
      for j = 1, numVals do
         local cItr = compFields[j]:get(cIndexer(idx))
         for i = 1, basis:numBasis() do
            assert_equal(vItr[i+(j-1)*basis:numBasis()], cItr[i], "Checking components")
         end
      end
   end
end

-- Run tests.
test()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
