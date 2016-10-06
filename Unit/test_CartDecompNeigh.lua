-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition neighbor calculations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local DecompRegionCalc = require "Lib.CartDecomp"
local CartDecompNeigh = require "Lib.CartDecompNeigh"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local Time = require "Lib.Time"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function hasVal(tbl, val)
   for _,v in ipairs(tbl) do
      if v == val then return true end
   end
   return false
end

function test_1()
   -- create decomposition and decompose a region
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
   decomp:decompose(Range.Range({1, 1}, {10, 10}))
   -- create neighbor calculator
   local decompNeigh = CartDecompNeigh(decomp)

   -- decomp (0,0)
   decompNeigh:calcFaceCommNeigh(0, 0)
   for i = 1, decomp:numSubDomains() do
      local neighData = decompNeigh:neighborData(i)
      assert_equal(0, #neighData, "Checking that there should be no neigbors")
   end

   -- decomp (0,1)
   decompNeigh:calcFaceCommNeigh(0, 1)
   local neighData = decompNeigh:neighborData(1)
   assert_equal(true, hasVal(neighData, 2), "Testing decomp 1->2")
   assert_equal(true, hasVal(neighData, 3), "Testing decomp 1->3")
   assert_equal(false, hasVal(neighData, 4), "Testing decomp 1->4")

   local neighData = decompNeigh:neighborData(2)
   assert_equal(false, hasVal(neighData, 1), "Testing decomp 2->1")
   assert_equal(true, hasVal(neighData, 4), "Testing decomp 2->4")
   assert_equal(false, hasVal(neighData, 3), "Testing decomp 2->3")

   local neighData = decompNeigh:neighborData(3)
   assert_equal(false, hasVal(neighData, 1), "Testing decomp 3->1")
   assert_equal(true, hasVal(neighData, 4), "Testing decomp 3->4")
   assert_equal(false, hasVal(neighData, 2), "Testing decomp 3->2")

   local neighData = decompNeigh:neighborData(4)
   assert_equal(false, hasVal(neighData, 3), "Testing decomp 4->3")
   assert_equal(false, hasVal(neighData, 2), "Testing decomp 4->2")
   assert_equal(false, hasVal(neighData, 1), "Testing decomp 4->1")

   -- decomp (1,1)
   decompNeigh:calcFaceCommNeigh(1, 1)
   local neighData = decompNeigh:neighborData(1)
   assert_equal(true, hasVal(neighData, 2), "Testing decomp")
   assert_equal(true, hasVal(neighData, 3), "Testing decomp")
   assert_equal(false, hasVal(neighData, 4), "Testing decomp")

   local neighData = decompNeigh:neighborData(2)
   assert_equal(true, hasVal(neighData, 1), "Testing decomp")
   assert_equal(true, hasVal(neighData, 4), "Testing decomp")
   assert_equal(false, hasVal(neighData, 3), "Testing decomp")

   local neighData = decompNeigh:neighborData(3)
   assert_equal(true, hasVal(neighData, 1), "Testing decomp")
   assert_equal(true, hasVal(neighData, 4), "Testing decomp")
   assert_equal(false, hasVal(neighData, 2), "Testing decomp")

   local neighData = decompNeigh:neighborData(4)
   assert_equal(true, hasVal(neighData, 3), "Testing decomp")
   assert_equal(true, hasVal(neighData, 2), "Testing decomp")
   assert_equal(false, hasVal(neighData, 1), "Testing decomp")
end

function test_2()
   -- create decomposition and decompose a region
   local decomp = DecompRegionCalc.CartProd { cuts = {3, 3} }
   decomp:decompose(Range.Range({1, 1}, {30, 30}))
   -- create neighbor calculator
   local decompNeigh = CartDecompNeigh(decomp)

   decompNeigh:calcFaceCommNeigh(1, 1)
   for i = 1, decomp:numSubDomains() do
      local nd = decompNeigh:neighborData(i)
      if i==1 or i==3 or i==9 or i==7  then
	 assert_equal(2, #nd, string.format("Testing size of neighbors for dom=%d", i))
      elseif i==2 or i==6 or i==8 or i==4 then
	 assert_equal(3, #nd, string.format("Testing size of neighbors for dom=%d", i))
      else
	 assert_equal(4, #nd, "Testing size of neighbors for dom=5")
      end
   end
end

function test_3()
   local tmStart = Time.clock()   
   -- create decomposition and decompose a region
   local decomp = DecompRegionCalc.CartProd { cuts = {20, 20, 20} }
   decomp:decompose(Range.Range({1, 1, 1}, {2000, 2000, 2000}))
   -- create neighbor calculator
   local decompNeigh = CartDecompNeigh(decomp)

   decompNeigh:calcFaceCommNeigh(1, 1)
   local tmCalc = Time.clock()-tmStart
   print(string.format("test_3() took %g secs for %d sub-domains", tmCalc, decomp:numSubDomains()))
end

-- Run tests
test_3()
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
