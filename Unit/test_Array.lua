-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Array = require "DataStruct.Array"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local arr = Array.Array( {10, 20}, Array.double)

   assert_equal(2, arr:rank(), "Testing array rank")
   local shape = arr:shape()
   assert_equal(10, shape[1], "Testing array shape")
   assert_equal(20, shape[2], "Testing array shape")
   assert_equal(1, arr.c, "Use count")
   assert_equal(200, arr:size(), "Total n")
   assert_equal(ffi.sizeof("double"), arr:elemSize(), "Size of element")

   local brr = arr:clone()
   assert_equal(2, brr:rank(), "Testing brray rank")
   shape = brr:shape()
   assert_equal(10, shape[1], "Testing brray shape")
   assert_equal(20, shape[2], "Testing brray shape")
   assert_equal(1, brr.c, "Use count")
   assert_equal(200, brr:size(), "Total n")
   assert_equal(ffi.sizeof("double"), brr:elemSize(), "Size of element")
end

function test_2()
   local arr = Array.Array( {10}, Array.double)
   -- fetch pointer to data
   local arr_d = ffi.cast("double *", arr.d)

   -- set values
   for i = 0, tonumber(arr.s[0])-1 do
      arr_d[i] = i+0.5
   end

   local brr = arr:clone()
   local brr_d = ffi.cast("double *", brr.d)
   -- check values
   for i = 0, tonumber(brr.s[0])-1 do
      assert_equal(i+0.5, brr_d[i], "Checking clone")
   end

   local crr = Array.Array(arr:shape(), Array.double)
   crr:copy(arr)
   crr_d = ffi.cast("double *", crr.d)
   for i = 0, tonumber(brr.s[0])-1 do
      assert_equal(arr_d[i], crr_d[i], "Checking copy " .. i)
   end
end

function test_3()
   local a = Array.Array( {10, 15, 20}, Array.float)
   
   -- Aquire pointer and store in a C-struct
   local ds = ffi.new(ffi.typeof("struct { GkylArray_t *d; } "), a:aquire())
   assert_equal(2, a.c, "Checking use-count")

   -- Original 'a' can go away at this point
   a = nil
   collectgarbage("collect") -- force 'a' to dissapear

   assert_equal(1, ds.d.c, "Checking use-count")
   assert_equal(3, ds.d:rank(), "Rank")

   -- Release pointer to prevent memory leak
   ds.d:release()
end

function test_4()
   local a = Array.Array( {10, 15, 20}, Array.float)
   
   -- Aquire pointer and store in a C-struct: in this test we attach a
   -- destructor to the created object in which we release() the
   -- aquire()-ed pointer
   local ds = ffi.gc(
      ffi.new(ffi.typeof("struct { GkylArray_t *d; } "), a:aquire()),
      function (p) p.d:release() end
   )
   assert_equal(2, a.c, "Checking use-count")

   -- Original 'a' can go away at this point
   a = nil
   collectgarbage("collect") -- force 'a' to dissapear

   assert_equal(1, ds.d.c, "Checking use-count")
   assert_equal(3, ds.d:rank(), "Rank")
end

function test_5()
   local a = Array.Array( {10},
      ffi.typeof("struct { double x, y, z;}")
   )
   assert_equal(3*ffi.sizeof("double"), a:elemSize(), "Checking size of custom type")

   local dptr = ffi.cast(ffi.typeof("struct { double x, y, z;} *"), a.d)
   dptr[0].x, dptr[0].y, dptr[0].z = 1.0, 2.0, 3.0
end

function test_6()
   local a = Array.Array( {10, 15, 20}, Array.int)

   a:reshape { 10*15, 20 }

   assert_equal(2, a:rank(), "Checking new rank")
   local shape = a:shape()
   assert_equal(10*15, shape[1], "Checking new shape")
   assert_equal(20, shape[2], "Checking new shape")

   a:reshape { 10, 3, 5, 4, 5 }
   assert_equal(5, a:rank(), "Checking new rank")
   shape = a:shape()
   assert_equal(10, shape[1], "Checking new shape")
   assert_equal(3, shape[2], "Checking new shape")
   assert_equal(5, shape[3], "Checking new shape")
   assert_equal(4, shape[4], "Checking new shape")
   assert_equal(5, shape[5], "Checking new shape")
end

-- Run tests
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
