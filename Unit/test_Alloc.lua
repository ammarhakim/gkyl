-- Gkyl ------------------------------------------------------------------------
--
-- Test for memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Alloc = require "Lib.Alloc"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local len = 100
   local da = Alloc.Double(len)
   da:fill(1.0)
   for i = 1, len do
      assert_equal(1.0, da[i], "Checking fill")
   end
   for i = 1, len do
      da[i] = (i+0.5)*0.1
   end
   for i = 1, len do
      assert_equal((i+0.5)*0.1, da[i], "Checking [] operator")
   end
end

function test_2()
   local len = 100
   local eulerAlloc = Alloc.createAllocator("struct {double rho, rhou, E;}")
   local eulerFld = eulerAlloc(len)

   for i = 1, len do
      eulerFld[i].rho = (i+0.5)*0.1
      eulerFld[i].rhou = (i+0.5)*0.2
      eulerFld[i].E = (i+0.5)*0.3
   end

   for i = 1, len do
      assert_equal((i+0.5)*0.1, eulerFld[i].rho, "Checking [] operator")
      assert_equal((i+0.5)*0.2, eulerFld[i].rhou, "Checking [] operator")
      assert_equal((i+0.5)*0.3, eulerFld[i].E, "Checking [] operator")
   end   
end

-- run tes
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
