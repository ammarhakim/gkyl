-- Gkyl ------------------------------------------------------------------------
--
-- Test for SQLITE3 library
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local sql = require "sqlite3"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1(conn)

   -- create table
   conn:exec[[ CREATE TABLE t(id TEXT, num REAL); ]]

   -- insert data
   for i = 1, 4 do
      conn:exec(string.format("insert into t values('%s', %d)", 'var'..i, 200*i))
   end

   -- festch data and check
   local t, nrow = conn:exec("SELECT * FROM t") -- Records are by column.
   local tId, tNum = t["id"], t["num"]
   for i = 1, nrow do
      assert_equal('var'..i, tId[i], "Checking row ID")
      assert_equal(200*i, tNum[i], "Checking row value")
   end
end

-- Run tests
local conn = sql.open("") -- temporary in-memory database
test_1(conn)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
