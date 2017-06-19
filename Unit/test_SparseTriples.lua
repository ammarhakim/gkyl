-- Gkyl ------------------------------------------------------------------------
--
-- Test for sparse triples
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local st = DataStruct.SparseTriples { 10, 20 }

   assert_equal(10, st:numRows(), "Testing number of rows")
   assert_equal(20, st:numColumns(), "Testing number of cols")

   -- add some triples
   st:append(1, 1, 34.5)
   st:append(1, 2, 44.5)

   assert_equal(2, st:size(), "Testing size")

   local t = st:g(1)
   assert_equal(1, t.i, "Testing triple")
   assert_equal(1, t.j, "Testing triple")
   assert_equal(34.5, t.val, "Testing triple")

   t = st:g(2)
   assert_equal(1, t.i, "Testing triple")
   assert_equal(2, t.j, "Testing triple")
   assert_equal(44.5, t.val, "Testing triple")   
end

-- run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
