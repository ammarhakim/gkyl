-- Gkyl ------------------------------------------------------------------------
--
-- Test for SQLITE3 library
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local sql = require "sqlite3"
local Time = require "Lib.Time"
local date = require "xsys.date"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1(conn)
   -- create table
   conn:exec[[ CREATE TABLE t(id TEXT, num REAL); ]]

   -- insert data
   for i = 1, 4 do
      conn:exec(string.format("insert into t values('%s', %d)", 'var'..i, 200*i))
   end

   -- fetch data and check
   local t, nrow = conn:exec("SELECT * FROM t") -- Records are by column.
   local tId, tNum = t["id"], t["num"]
   for i = 1, nrow do
      assert_equal('var'..i, tId[i], "Checking row ID")
      assert_equal(200*i, tNum[i], "Checking row value")
   end
end

function test_2(conn)
   conn:exec [[ 
    create table RegressionData (
      tstamp text,
      name text,
      GKYL_EXEC text,
      GKYL_HG_CHANGESET text,
      GKYL_BUILD_DATE text,
      status integer,
      runtime real
    )
   ]]

   -- stored procedure to insert data into table
   local stmt = conn:prepare [[
      insert into RegressionData values (
        ?, ?, ?, ?, ?, ?, ?
      )
   ]]
   local function insertData(tm, nm, status, runTm)
      stmt:reset():bind(
	 tm,
	 nm,
	 GKYL_EXEC,
	 GKYL_GIT_CHANGESET,
	 GKYL_BUILD_DATE,
	 status == true and 1 or 0,
	 runTm):step()
   end
   insertData(date(false):fmt(), "rt-two-stream-p1.lua", true, 1.5)
   insertData(date(false):fmt(), "rt-two-stream-p2.lua", false, 0.0)
   insertData(date(false):fmt(), "rt-two-stream-p1.lua", true, 10.5)
   insertData(date(false):fmt(), "rt-two-stream-p1.lua", true, 3.5)
   insertData(date(false):fmt(), "rt-two-stream-p2.lua", true, 3.1)

   local t, nrow
      = conn:exec("select runtime from RegressionData where name=='rt-two-stream-p1.lua'")

   assert_equal(3, nrow, "Checking number of rows in table for rt-two-stream-p1.lua")
   assert_equal(1.5, t["runtime"][1], "Checking runtime column")
   assert_equal(10.5, t["runtime"][2], "Checking runtime column")
   assert_equal(3.5, t["runtime"][3], "Checking runtime column")   
end

-- Run tests
local conn = sql.open("") -- temporary in-memory database
test_1(conn)
test_2(conn)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
