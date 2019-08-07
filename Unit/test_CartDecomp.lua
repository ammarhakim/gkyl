-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition of Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local CartDecomp = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local decomp = CartDecomp.CartProd { cuts = {2, 3}, __serTesting = true }
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
   local decomp = CartDecomp.CartProd { cuts = {2, 3, 4}, __serTesting = true }
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
   local decomp = CartDecomp.CartProd { cuts = {1, 1}, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))

   local perIds = decomposedRgn:boundarySubDomainIds(1)
   assert_equal(1, #perIds, "Checking number of skeleton sub-domains in direction 1")
   -- check individual pairs
   assert_equal(1, perIds[1].lower, "Checking lower")
   assert_equal(1, perIds[1].upper, "Checking upper")

   perIds = decomposedRgn:boundarySubDomainIds(2)
   assert_equal(1, #perIds, "Checking number of skeleton sub-domains in direction 2")
   -- check individual pairs
   assert_equal(1, perIds[1].lower, "Checking lower")
   assert_equal(1, perIds[1].upper, "Checking upper")   
end

function test_4()
   local decomp = CartDecomp.CartProd { cuts = {2, 3}, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))

   local perIds = decomposedRgn:boundarySubDomainIds(1)
   assert_equal(3, #perIds, "Checking number of skeleton sub-domains in direction 1")
   -- check individual pairs
   assert_equal(1, perIds[1].lower, "Checking lower")
   assert_equal(2, perIds[1].upper, "Checking upper")

   assert_equal(3, perIds[2].lower, "Checking lower")
   assert_equal(4, perIds[2].upper, "Checking upper")

   assert_equal(5, perIds[3].lower, "Checking lower")
   assert_equal(6, perIds[3].upper, "Checking upper")   


   local perIds = decomposedRgn:boundarySubDomainIds(2)
   assert_equal(2, #perIds, "Checking number of skeleton sub-domains in direction 2")
   -- check individual pairs
   assert_equal(1, perIds[1].lower, "Checking lower")
   assert_equal(5, perIds[1].upper, "Checking upper")

   assert_equal(2, perIds[2].lower, "Checking lower")
   assert_equal(6, perIds[2].upper, "Checking upper")
end

function showRange(msg, rng)
   print(msg)
   print(string.format("  (%d, %d) - (%d, %d). Vol = %d",
		       rng:lower(1), rng:lower(2), rng:upper(1), rng:upper(2), rng:volume()))
end

function test_5()
   local decomp = CartDecomp.CartProd { cuts = {2, 3}, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {20, 30}))

   for d = 1, decomp:ndim() do
      -- fetch skeleton in direction d
      local perIds = decomposedRgn:boundarySubDomainIds(d)
      for p = 1, #perIds do
	 local rlo, rup = decomposedRgn:subDomain(perIds[p].lower), decomposedRgn:subDomain(perIds[p].upper)
	 --showRange(string.format("Dir %d. Lower-side box", d), rlo)
	 --showRange(string.format("Dir %d. Upper-side box", d), rup)
      end
   end

end

function test_6()
   ---
   local cells = {10, 10}
   local cuts = CartDecomp.makeCuts(#cells, 4, cells)
   
   assert_equal(2, cuts[1], "Checking cuts")
   assert_equal(2, cuts[2], "Checking cuts")

   ----
   local cells = {10, 20}
   local cuts = CartDecomp.makeCuts(#cells, 6, cells)
   
   assert_equal(2, cuts[1], "Checking cuts")
   assert_equal(3, cuts[2], "Checking cuts")

   --
   local cells = {20, 10}
   local cuts = CartDecomp.makeCuts(#cells, 6, cells)
   
   assert_equal(3, cuts[1], "Checking cuts")
   assert_equal(2, cuts[2], "Checking cuts")

   --
   local cells = {20, 10, 30}
   local cuts = CartDecomp.makeCuts(#cells, 6, cells)
   
   assert_equal(2, cuts[1], "Checking cuts")
   assert_equal(1, cuts[2], "Checking cuts")
   assert_equal(3, cuts[3], "Checking cuts")

   --
   local cells = {20, 10, 30}
   local cuts = CartDecomp.makeCuts(#cells, 1234567, cells)
   
   assert_equal(127, cuts[1], "Checking cuts")
   assert_equal(1, cuts[2], "Checking cuts")
   assert_equal(9721, cuts[3], "Checking cuts")

   --
   local cells = {20, 20, 20}
   local cuts = CartDecomp.makeCuts(#cells, 1000, cells)
   
   assert_equal(10, cuts[1], "Checking cuts")
   assert_equal(10, cuts[2], "Checking cuts")
   assert_equal(10, cuts[3], "Checking cuts")

   --
   local cells = {20, 20, 20}
   local cuts = CartDecomp.makeCuts(#cells, 100^3, cells)
   
   assert_equal(100, cuts[1], "Checking cuts")
   assert_equal(100, cuts[2], "Checking cuts")
   assert_equal(100, cuts[3], "Checking cuts")   
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
