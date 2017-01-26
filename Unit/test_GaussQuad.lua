-- Gkyl
-- ------------------------------------------------------------------------
--
-- Test for Guass quadrature rules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local GaussQuadRules = require "Lib.GaussQuadRules"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function sumArray(n, arr)
   local sum = 0.0
   for i = 1, n do
      sum = sum + arr[i]
   end
   return sum
end


function integrate(n, mu, w, func)
   local sum = 0.0
   for i = 1, n do
      sum = sum + w[i]*func(mu[i])
   end
   return sum
end

function test_1()
   for i = 1, #GaussQuadRules.ordinates do
      assert_equal(2.0, sumArray(i, GaussQuadRules.weights[i], "Testing to make sure weights sum to 2. N = " .. i))
   end

   local muList, wList = GaussQuadRules.ordinates, GaussQuadRules.weights

   -- test integrations: all of these are exact
   
   local N=1   
   local r = integrate(N, muList[N], wList[N], function (x) return 1 end)
   assert_equal(r, 2, "Testing integral of 1")

   N=2
   r = integrate(N, muList[N], wList[N], function (x) return x^2 end)
   assert_equal(r, 2/3, "Testing integral of x^2")

   N=3
   r = integrate(N, muList[N], wList[N], function (x) return x^4+x^2 end)
   assert_equal(r, 16/15, "Testing integral of x^4+x^2")

   N=4
   r = integrate(N, muList[N], wList[N], function (x) return x^6+x^4+x^2 end)
   assert_equal(r, 142/105, "Testing integral of x^6+x^4+x^2")

   N=5
   r = integrate(N, muList[N], wList[N], function (x) return x^8+x^6+x^4+x^2 end)
   assert_equal(r, 496/315, "Testing integral of x^8+x^6+x^4+x^2")

   N=6
   r = integrate(N, muList[N], wList[N], function (x) return x^10+x^8+x^6+x^4+x^2 end)
   assert_equal(r, 6086/3465, "Testing integral of x^10+x^8+x^6+x^4+x^2")

   N=7
   r = integrate(N, muList[N], wList[N], function (x) return x^12+x^10+x^8+x^6+x^4+x^2 end)
   assert_equal(r, 86048/45045, "Testing integral of x^12+x^10+x^8+x^6+x^4+x^2")

   N=8
   r = integrate(N, muList[N], wList[N], function (x) return x^14+x^12+x^10+x^8+x^6+x^4+x^2 end)
   assert_equal(r, 92054/45045, "Testing integral of x^14+x^12+x^10+x^8+x^6+x^4+x^2")
   
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
