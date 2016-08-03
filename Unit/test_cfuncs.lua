-- Gkyl ------------------------------------------------------------------------
--
-- Test for LuaJIT/C bridge
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local xsys = require "xsys"
local Grid = require "Grid"
local DataStruct = require "DataStruct"

local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  typedef struct { double x, y, z; } loc_t;
  double calcSum(int n, double *v);
  double addValues(loc_t *v);
  void setValues(int n, int ix, double *v);
]]

function test_1()
   local v = ffi.new(typeof("double[?]"), 10)
   for i = 1, 10 do
      v[i-1] = i
   end
   local sum = ffi.C.calcSum(10, v)
   assert_equal(55, sum, "Checking if external call to sum worked")
end

function test_2()
   local v = ffi.new(typeof("loc_t"))
   v.x, v.y, v.z = 1.0, 2.0, 3.0
   local sum = ffi.C.addValues(v)
   assert_equal(6, sum, "Checking if external call to sum worked")
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {10},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
   }

   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      ffi.C.setValues(3, i, fitr._cdata)
   end

   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      assert_equal(i+1, fitr[1], "Checking field value")
      assert_equal(i+2, fitr[2], "Checking field value")
      assert_equal(i+3, fitr[3], "Checking field value")
   end
end

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
