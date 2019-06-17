-- Gkyl ------------------------------------------------------------------------
--
-- Test for prime factor algorithm
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Lin = require "Lib.Linalg"
local PrimeFactor = require "Lib.PrimeFactor"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   -- 1 = {1}
   local v = PrimeFactor.largest(1)
   assert_equal(1, v, "Checking largest prime-factor")

   local av = PrimeFactor.all(1)
   assert_equal(1, av[1], "Checking prime-factor")

   -- 1000 = {2, 2, 2, 5, 5}
   local v = PrimeFactor.largest(1000)
   assert_equal(5, v, "Checking largest prime-factor")

   local av = PrimeFactor.all(1000)
   assert_equal(2, av[1], "Checking prime-factor")
   assert_equal(2, av[2], "Checking prime-factor")
   assert_equal(2, av[3], "Checking prime-factor")
   assert_equal(5, av[4], "Checking prime-factor")
   assert_equal(5, av[5], "Checking prime-factor")
   
   -- 1234567 = {127, 9721}
   local v = PrimeFactor.largest(1234567)
   assert_equal(9721, v, "Checking largest prime-factor")

   local av = PrimeFactor.all(1234567)
   assert_equal(127, av[1], "Checking prime-factor")
   assert_equal(9721, av[2], "Checking prime-factor")

   -- 12345678901 = {857, 14405693}
   local v = PrimeFactor.largest(12345678901)
   assert_equal(14405693, v, "Checking largest prime-factor")

   local av = PrimeFactor.all(12345678901)
   assert_equal(857, av[1], "Checking prime-factor")
   assert_equal(14405693, av[2], "Checking prime-factor")   
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
