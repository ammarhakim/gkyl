-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition of Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 3}, __serTesting = true }
   assert_equal(2, decomp:ndim(), "Checking ndim")
   
   assert_equal(2, decomp:cuts(1), "Checking cuts")
   assert_equal(3, decomp:cuts(2), "Checking cuts")

   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))
   assert_equal(6, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- fetch domains and do sanity checks
   local v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):volume()
   end
   assert_equal(100, v, "Checking volume of decomposedRgn")

   v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(3*10, v, "Checking total number of X cells")

   v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):shape(2)
   end
   assert_equal(2*10, v, "Checking total number of Y cells")
   
end

function test_2()
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 3, 4}, __serTesting = true }
   assert_equal(3, decomp:ndim(), "Checking ndim")
   
   assert_equal(2, decomp:cuts(1), "Checking cuts")
   assert_equal(3, decomp:cuts(2), "Checking cuts")
   assert_equal(4, decomp:cuts(3), "Checking cuts")

   local decomposedRgn = decomp:decompose(Range.Range({1, 2, 3}, {10, 11, 12}))
   assert_equal(24, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- fetch domains and do sanity checks
   local v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):volume()
   end
   assert_equal(1000, v, "Checking volume of decomposedRgn")
   
   v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(12*10, v, "Checking total number of X cells")

   v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):shape(2)
   end
   assert_equal(8*10, v, "Checking total number of Y cells")
   
end

function test_3()
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 3}, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))

   local perIds = decomposedRgn:boundarySubDomainIds(1)
   assert_equal(3, #perIds, "Checking number of skeleton sub-domains in direction 1")

   local perIds = decomposedRgn:boundarySubDomainIds(2)
   assert_equal(2, #perIds, "Checking number of skeleton sub-domains in direction 2")
end

-- Run tests
test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
