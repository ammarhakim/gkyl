local Grid = require "Grid"
local DataStruct = require "DataStruct"
local StairSteppedBc = require "Updater.StairSteppedBc"
local BoundaryCondition = require "Updater.BoundaryCondition"

local grid = Grid.RectCart {
   lower = {0.},
   upper = {1.},
   cells = {10},
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
      local fieldItr = field:get(idxr(i))
      local inOutItr = inOut:get(inOutIdxr(i))
      print("i", i, "field", fieldItr[1], "inOut", inOutItr[1])
   end
end

local init = function()
   for i = localRange:lower(1), localRange:upper(1) do
      local fieldItr = field:get(idxr(i))
      local inOutItr = inOut:get(inOutIdxr(i))
      fieldItr[1] = 1.
      inOutItr[1] = 1
      if (i>3 and i<7) then
         inOutItr[1] = -1
      end
   end
end

local test = function(bcList)
   local ssBc = StairSteppedBc {
      onGrid = grid,
      inOut = inOut,
      boundaryConditions = bcList,
      dir = 1,
   }
   ssBc:_advance(0, 0, {field}, {field})
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
