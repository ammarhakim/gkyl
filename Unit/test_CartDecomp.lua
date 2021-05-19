-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition of Cartesian grids
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local CartDecomp = require "Lib.CartDecomp"
local Lin        = require "Lib.Linalg"
local Range      = require "Lib.Range"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"


local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1()
   local decomp = CartDecomp.CartProd { cuts = {2, 3}, __serTesting = true }
   assert_equal(2, decomp:ndim(), "Checking ndim")
   
   assert_equal(2, decomp:cuts(1), "Checking X cuts")
   assert_equal(3, decomp:cuts(2), "Checking Y cuts")

   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))
   assert_equal(6, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- Fetch domains and do sanity checks.
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
   
   assert_equal(2, decomp:cuts(1), "Checking X cuts")
   assert_equal(3, decomp:cuts(2), "Checking Y cuts")
   assert_equal(4, decomp:cuts(3), "Checking Z cuts")

   local decomposedRgn = decomp:decompose(Range.Range({1, 2, 3}, {10, 11, 12}))
   assert_equal(24, decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   -- Fetch domains and do sanity checks.
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
   local decomp        = CartDecomp.CartProd { cuts = {1, 1}, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))

   local perIds = decomposedRgn:boundarySubDomainIds(1)
   assert_equal(1, #perIds, "Checking number of skeleton sub-domains in direction 1")
   -- Check individual pairs.
   assert_equal(1, perIds[1].lower, "Checking X lower")
   assert_equal(1, perIds[1].upper, "Checking X upper")

   perIds = decomposedRgn:boundarySubDomainIds(2)
   assert_equal(1, #perIds, "Checking number of skeleton sub-domains in direction 2")
   -- Check individual pairs.
   assert_equal(1, perIds[1].lower, "Checking Y lower")
   assert_equal(1, perIds[1].upper, "Checking Y upper")   
end

function test_4()
--   local decomp        = CartDecomp.CartProd { cuts = {2, 1}, __serTesting = true }
   local decomp        = CartDecomp.CartProd { cuts = {2, 1} }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, {10, 10}))

   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 10},
      decomposition = decomp,
      periodicDirs = {2},
   }
   local field = DataStruct.Field {
      onGrid        = grid,
      numComponents = 1,
      ghost         = {1, 1},
      syncCorners   = true,
   }

--   local perIds = decomposedRgn:boundarySubDomainIds(1)
--   assert_equal(3, #perIds, "Checking number of skeleton sub-domains in direction 1")
--   -- Check individual pairs.
--   assert_equal(1, perIds[1].lower, "Checking X lower")
--   assert_equal(2, perIds[1].upper, "Checking X upper")
--
--   assert_equal(3, perIds[2].lower, "Checking X lower")
--   assert_equal(4, perIds[2].upper, "Checking X upper")
--
--   assert_equal(5, perIds[3].lower, "Checking X lower")
--   assert_equal(6, perIds[3].upper, "Checking X upper")   
--
--
--   local perIds = decomposedRgn:boundarySubDomainIds(2)
--   assert_equal(2, #perIds, "Checking number of skeleton sub-domains in direction 2")
--   -- Check individual pairs.
--   assert_equal(1, perIds[1].lower, "Checking Y lower")
--   assert_equal(5, perIds[1].upper, "Checking Y upper")
--
--   assert_equal(2, perIds[2].lower, "Checking Y lower")
--   assert_equal(6, perIds[2].upper, "Checking Y upper")
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
      -- Fetch skeleton in direction d.
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

function test_7()
   -- Test the :childDecomp method with a 1D decomp.
   local cuts1D  = {3}
   local cells1D = {10}
   local decomp  = CartDecomp.CartProd { cuts = cuts1D, __serTesting = true }

   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({1})

   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }

   assert_equal(1, childDecomp:ndim(), "Checking ndim")
   
   assert_equal(cuts1D[1], childDecomp:cuts(1), "Checking cuts")

   local decomposedRgn = childDecomp:decompose(Range.Range({1}, cells1D))
   assert_equal(cuts1D[1], decomposedRgn:numSubDomains(), "Checking number of sub-domains")

   local v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):volume()
   end
   assert_equal(cells1D[1], v, "Checking volume of decomposedRgn")

   v = 0
   for i = 1, decomposedRgn:numSubDomains() do
      v = v + decomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells1D[1], v, "Checking total number of X cells")

   -- Test the :childDecomp method with a 2D decomp.
   local cuts2D  = {2, 3}
   local cells2D = {10, 48}
   local decomp = CartDecomp.CartProd { cuts = cuts2D, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1}, cells2D))

   -- Sub-decomp along first dimension.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({1})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(1, childDecomp:ndim(), "Checking x-child ndim")
   assert_equal(cuts2D[1], childDecomp:cuts(1), "Checking x-child cuts")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1}, {cells2D[1]}))
   assert_equal(cuts2D[1], childDecomposedRgn:numSubDomains(), "Checking number of x-child sub-domains")
   local v, w = 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells2D[1], v, "Checking volume of childDecomposedRgn along 1st dimension")
   assert_equal(cells2D[1], w, "Checking total number of child X cells")

   -- Sub-decomp along second dimension.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({2})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(1, childDecomp:ndim(), "Checking y-child ndim")
   assert_equal(cuts2D[2], childDecomp:cuts(1), "Checking y-child cuts")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1}, {cells2D[2]}))
   assert_equal(cuts2D[2], childDecomposedRgn:numSubDomains(), "Checking number of y-child sub-domains")
   local v, w = 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells2D[2], v, "Checking volume of childDecomposedRgn along 2nd dimension")
   assert_equal(cells2D[2], w, "Checking total number of child Y cells")


   -- Test the :childDecomp method with a 3D decomp.
   local cuts3D  = {2, 8, 5}
   local cells3D = {10, 48, 15}
   local decomp = CartDecomp.CartProd { cuts = cuts3D, __serTesting = true }
   local decomposedRgn = decomp:decompose(Range.Range({1, 1, 1}, cells3D))

   -- Sub-decomp along 1st dimension.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({1})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(1, childDecomp:ndim(), "Checking x-child ndim of 3D decomp")
   assert_equal(cuts3D[1], childDecomp:cuts(1), "Checking x-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1}, {cells3D[1]}))
   assert_equal(cuts3D[1], childDecomposedRgn:numSubDomains(), "Checking number of x-child sub-domains of 3D decomp")
   local v, w = 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells3D[1], v, "Checking volume of childDecomposedRgn along 1st dimension of 3D decomp")
   assert_equal(cells3D[1], w, "Checking total number of child X cells of 3D decomp")

   -- Sub-decomp along 2nd dimension.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({2})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(1, childDecomp:ndim(), "Checking y-child ndim of 3D decomp")
   assert_equal(cuts3D[2], childDecomp:cuts(1), "Checking y-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1}, {cells3D[2]}))
   assert_equal(cuts3D[2], childDecomposedRgn:numSubDomains(), "Checking number of y-child sub-domains of 3D decomp")
   local v, w = 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells3D[2], v, "Checking volume of childDecomposedRgn along 2nd dimension of 3D decomp")
   assert_equal(cells3D[2], w, "Checking total number of child Y cells of 3D decomp")

   -- Sub-decomp along 3rd dimension.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({3})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(1, childDecomp:ndim(), "Checking z-child ndim of 3D decomp")
   assert_equal(cuts3D[3], childDecomp:cuts(1), "Checking z-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1}, {cells3D[3]}))
   assert_equal(cuts3D[3], childDecomposedRgn:numSubDomains(), "Checking number of z-child sub-domains of 3D decomp")
   local v, w = 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
   end
   assert_equal(cells3D[3], v, "Checking volume of childDecomposedRgn along 3rd dimension of 3D decomp")
   assert_equal(cells3D[3], w, "Checking total number of child Z cells of 3D decomp")

   -- Sub-decomp along 1st and 2nd dimensions.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({1,2})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(2, childDecomp:ndim(), "Checking xy-child ndim of 3D decomp")
   assert_equal(cuts3D[1], childDecomp:cuts(1), "Checking xy-child cuts of 3D decomp")
   assert_equal(cuts3D[2], childDecomp:cuts(2), "Checking xy-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1, 1}, {cells3D[1], cells3D[2]}))
   assert_equal(cuts3D[1]*cuts3D[2], childDecomposedRgn:numSubDomains(), "Checking number of xy-child sub-domains of 3D decomp")
   local v, w, x = 0, 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
      x = x + childDecomposedRgn:subDomain(i):shape(2)
   end
   assert_equal(cells3D[1]*cells3D[2], v, "Checking volume of xy childDecomposedRgn of 3D decomp")
   assert_equal(cells3D[1]*cuts3D[2], w, "Checking total number of xy-child X cells of 3D decomp")
   assert_equal(cells3D[2]*cuts3D[1], x, "Checking total number of xy-child Y cells of 3D decomp")

   -- Sub-decomp along 2nd and 3rd dimensions.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({2,3})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(2, childDecomp:ndim(), "Checking yz-child ndim of 3D decomp")
   assert_equal(cuts3D[2], childDecomp:cuts(1), "Checking yz-child cuts of 3D decomp")
   assert_equal(cuts3D[3], childDecomp:cuts(2), "Checking yz-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1, 1}, {cells3D[2], cells3D[3]}))
   assert_equal(cuts3D[2]*cuts3D[3], childDecomposedRgn:numSubDomains(), "Checking number of yz-child sub-domains of 3D decomp")
   local v, w, x = 0, 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
      x = x + childDecomposedRgn:subDomain(i):shape(2)
   end
   assert_equal(cells3D[2]*cells3D[3], v, "Checking volume of yz childDecomposedRgn of 3D decomp")
   assert_equal(cells3D[2]*cuts3D[3], w, "Checking total number of yz-child X cells of 3D decomp")
   assert_equal(cells3D[3]*cuts3D[2], x, "Checking total number of yz-child Y cells of 3D decomp")

   -- Sub-decomp along 1st and 3rd dimensions.
   local childComm, childWriteRank, childCuts, childIsShared = decomp:childDecomp({1,3})
   local childDecomp = CartDecomp.CartProd {
      comm      = childComm,
      writeRank = childWriteRank,
      cuts      = childCuts,
      useShared = childIsShared,
      __serTesting = true,
   }
   assert_equal(2, childDecomp:ndim(), "Checking xz-child ndim of 3D decomp")
   assert_equal(cuts3D[1], childDecomp:cuts(1), "Checking xz-child cuts of 3D decomp")
   assert_equal(cuts3D[3], childDecomp:cuts(2), "Checking xz-child cuts of 3D decomp")
   local childDecomposedRgn = childDecomp:decompose(Range.Range({1, 1}, {cells3D[1], cells3D[3]}))
   assert_equal(cuts3D[1]*cuts3D[3], childDecomposedRgn:numSubDomains(), "Checking number of xz-child sub-domains of 3D decomp")
   local v, w, x = 0, 0, 0
   for i = 1, childDecomposedRgn:numSubDomains() do
      v = v + childDecomposedRgn:subDomain(i):volume()
      w = w + childDecomposedRgn:subDomain(i):shape(1)
      x = x + childDecomposedRgn:subDomain(i):shape(2)
   end
   assert_equal(cells3D[1]*cells3D[3], v, "Checking volume of xz childDecomposedRgn of 3D decomp")
   assert_equal(cells3D[1]*cuts3D[3], w, "Checking total number of xz-child X cells of 3D decomp")
   assert_equal(cells3D[3]*cuts3D[1], x, "Checking total number of xz-child Y cells of 3D decomp")

end


-- Run tests
--test_1()
--test_2()
--test_3()
test_4()
--test_5()
--test_6()
--test_7()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
