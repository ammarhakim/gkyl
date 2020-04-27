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
local Grid = require "Grid"
local math = require("sci.math").generic
local Lin = require "Lib.Linalg"

local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  typedef struct { double x, y, z; } loc_t;
  double calcSum(int n, double *v);
  double addValues(loc_t *v);
  void setValues(int n, int ix, double *v);

  typedef struct Adder Adder;
  void *new_Adder(int n);
  int add_Adder(Adder* a, int x);
  int sub_Adder(Adder* a, int x);
  int setFuncPointer_Adder(Adder *a, int (*v)(int));
  int incr_Adder(Adder *a, int x);

  void setMetricFuncPointer_Adder(Adder *a, void (*gfunc)(double *xc, double *g));
  void printg_Adder(Adder *a, double r, double t);

  void *new_Dim_1();
  int getDim_Dim_1(void *d);
]]

ffi.cdef [[
  typedef struct Particle_type *Particle;
  typedef struct { Particle *p; } ptcl_t;
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

function test_4()
   local adder = ffi.C.new_Adder(10)
   assert_equal(14, ffi.C.add_Adder(adder, 4), "Testing value returned")
   assert_equal(6, ffi.C.sub_Adder(adder, 4), "Testing value returned")

   local function val(x) return 3*x-1 end
   
   ffi.C.setFuncPointer_Adder(adder, val)
   assert_equal(39, ffi.C.incr_Adder(adder, 10), "Func pointer test")
end

function test_5()
   local grid = Grid.MappedCart {
      lower = {1.0, 0.0},
      upper = {2.0, 2*math.pi},
      cells = {64, 64},
      -- map computational to physical space
      mapc2p = function(xc)
	 local r, theta = xc[1], xc[2]
	 return r*math.cos(theta), r*math.sin(theta)
      end
   }

   local adder = ffi.C.new_Adder(10)

   local myXc, myG = Lin.Vec(2), Lin.Vec(3)
   local function gfunc(xc, g)
      myXc[1], myXc[2] = xc[1], xc[2]
      grid:calcMetric(myXc, myG)
      g[1], g[2], g[3] = myG[1], myG[2], myG[3]
   end

   ffi.C.setMetricFuncPointer_Adder(adder, gfunc)

   ffi.C.printg_Adder(adder, 1.5, math.pi)
   
end

function test_6()
   local dim1 = ffi.C.new_Dim_1()
   assert_equal(1, ffi.C.getDim_Dim_1(dim1), "Checking if dim1 object is correct")
end

test_1()
test_2()
test_3()
test_4()
test_5()
test_6()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
