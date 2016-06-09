-- Gkyl ------------------------------------------------------------------------
--
-- Test for linear algebra objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local v = Lin.vec(3)
   assert_equal(3, #v, "Checking length of vector")
   -- set values
   for i= 1, #v do
      v[i] = (i+0.5)*0.1
   end
   -- test them
   for i= 1, #v do
      assert_equal((i+0.5)*0.1, v[i], "Checking vector content")
   end

   local vcopy = v:copy()
   -- test copy
   for i = 1, #vcopy do
      assert_equal((i+0.5)*0.1, vcopy[i], "Checking vector copy")
   end
end

function test_2()
   local eulerVec = Lin.new_vec_ct(ffi.typeof("struct {double rho, rhou, E;}"))
   local v = eulerVec(4)
   assert_equal(4, #v, "Checking length of vector")
   -- set values
   for i = 1, #v do
      v[i].rho = i+0.0
      v[i].rhou = i+1.0
      v[i].E = i+2.0 
   end
   -- test them
   for i = 1, #v do
      assert_equal(i+0.0, v[i].rho, "Checking vector of struct contents")
      assert_equal(i+1.0, v[i].rhou, "Checking vector of struct contents")
      assert_equal(i+2.0, v[i].E, "Checking vector of struct contents")
   end

   local vcopy = v:copy()
   -- test copy
   for i = 1, #vcopy do
      assert_equal(i+0.0, vcopy[i].rho, "Checking vector of struct contents copy")
      assert_equal(i+1.0, vcopy[i].rhou, "Checking vector of struct contents copy")
      assert_equal(i+2.0, vcopy[i].E, "Checking vector of struct contents copy")
   end
end

-- Run tests
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
