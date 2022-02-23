local Grid = require "Grid"
local DataStruct = require "DataStruct"
local StairSteppedBc = require "Updater.StairSteppedBc"
local BoundaryCondition = require "Updater.BoundaryCondition"

local grid = Grid.RectCart {
   lower = {0., 0.},
   upper = {1., 1.},
   cells = {10, 11},
}
local field = DataStruct.Field {
   onGrid = grid,
   numComponents = 1,
   ghost = {1, 1},
}
local inOut = DataStruct.Field {
   onGrid = grid,
   numComponents = 1,
   ghost = {1, 1},
}
local localRange = field:localRange()
local idxr = field:indexer()
local inOutIdxr = inOut:indexer()

local show = function()
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
         local fieldItr = field:get(idxr(i, j))
         local inOutItr = inOut:get(inOutIdxr(i, j))
         print("i", i, "j", j, "field", fieldItr[1], "inOut", inOutItr[1])
      end
   end
end

local init = function()
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
         local fieldItr = field:get(idxr(i, j))
         local inOutItr = inOut:get(inOutIdxr(i, j))
         fieldItr[1] = 1.
         inOutItr[1] = 1
         if (i>3 and i<7 and j>4 and j <8) then
            inOutItr[1] = -1
         end
      end
   end
end

local test = function(bcList)
   local ssBcX = StairSteppedBc {
      onGrid = grid,
      inOut = inOut,
      boundaryConditions = bcList,
      dir = 1,
   }
   local ssBcY = StairSteppedBc {
      onGrid = grid,
      inOut = inOut,
      boundaryConditions = bcList,
      dir = 2,
   }
   ssBcX:_advance(0, 0, {field}, {field})
   ssBcY:_advance(0, 0, {field}, {field})
end

init()
-- print("before ssBnd")
-- show()
field:write("field_before.bp", 0., 0, false)
inOut:write("inOut.bp", 0., 0, false)
test({
   BoundaryCondition.Const{components={1}, values={0.5}}
})
-- print("after ssBnd")
-- show()
field:write("field_after.bp", 0., 0, false)
