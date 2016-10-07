-- Gkyl ------------------------------------------------------------------------
--
-- Test for memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

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

-- run test
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
