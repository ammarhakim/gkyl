-- Gkyl ------------------------------------------------------------------------
--
-- Query code to get information from DB created by runregression code
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local sql = require "sqlite3"

-- load configuration file created by regression testing system
local f = loadfile("runregression.config.lua")
if not f then
   log("Regression tests not configured! Run config command first\n")
   os.exit(1)
end
configVals = f()

-- open connection
local sqlConn = sql.open(string.format("%s/regressiondb", configVals['results_dir']), "ro")

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("queryrdb")
   :require_command(true)
   :description [[
Query database created by runregression system to extract information
about various tests.
]]

-- read metadata table 
local function read_metatable()
   local t, nrow = sqlConn:exec("select * from RegressionMeta")
   local dbMeta = {}
   for i = 1, nrow do
      dbMeta[i] = {
	 guid = t['guid'][i],
	 tstamp = t['tstamp'][i],
	 changeset = t['GKYL_HG_CHANGESET'][i],
	 builddate = t['GKYL_BUILD_DATE'][i],
	 ntotal = t['ntotal'][i],
	 npass = t['npass'][i],
	 nfail = t['nfail'][i]
      }
   end
   return dbMeta
end

-- get all tests with specified GUID
local function read_tests_with_id(guid)
   local t, nrow = sqlConn:exec(
      string.format("select * from RegressionData where guid=='%s'", guid))
   local sstr = { [-2] = "create", [-1] = "skip", [0] = "fail", [1] = "pass" }

   local maxNm = 0
   
   dbData = {}
   for i = 1, nrow do
      dbData[i] = {
	 name = t['name'][i],
	 status = sstr[tonumber(t['status'][i])],
	 runtime = t['runtime'][i],
	 runlog = t['runlog'][i],
      }
      maxNm = math.max(maxNm, string.len(t['name'][i]))
   end
   return dbData, maxNm
end

-- function to handle summary ("s") command
local function summary_action(args, name)
   local dbMeta = read_metatable()
   local nrow = #dbMeta

   local fmt = "%-4s: %-20s %-30s %-5s %-5s %-5s"
   print(string.format(fmt, "ID", "Time-Stamp", "Changeset", "Total", "Pass", "Fail"))
   for i,d in pairs(dbMeta) do
      print(string.format(fmt, nrow-i+1, d.tstamp, d.changeset, tonumber(d.ntotal), tonumber(d.npass), tonumber(d.nfail)))
   end
end

-- function to handle query ("q") command
local function query_action(args, name)
   local dbMeta = read_metatable()
   local idx = #dbMeta-args.id+1 -- meta-list is printed in reverse order
   local dbData, maxNm = read_tests_with_id(dbMeta[idx].guid)

   if tonumber(args.test) > 0 then
      local tidx = math.min(args.test, #dbData)
      print(string.format("=== "))
      print(string.format("=== Log for test %s ==", dbData[tidx].name))
      print(string.format("=== ")) 
      print(dbData[tidx].runlog)
   else
      local fmt = "%-4s: %-" .. maxNm+2 .. "s %-7s %-4s"
      print(string.format(fmt, "ID", "Name", "Status", "Run-Time"))
      local fmt1 = "%-4s: %-" .. maxNm+2 .. "s %-7s %.4g"
      for i,d in pairs(dbData) do
	 print(string.format(fmt1, i, d.name, d.status, d.runtime))
      end
   end
end

-- function to handle query ("q") command
local function history_action(args, name)
   if args.regression == nil then return end

   local function trimname(nm)
      local s,e = string.find(nm, "./")
      if s and s==1 then return string.sub(nm, e+1) end
      return nm
   end

   local tNm = trimname(args.regression)
   local dat, nrow = sqlConn:exec(
      string.format("select * from RegressionData where name=='%s'", tNm)
   )

   local fmt = "%-20s %-30s %-4s"
   print(string.format(fmt, "Time-Stamp", "Changeset", "Run-Time"))
   local fmt1 = "%-20s %-30s %.4g"
   for i = 1,nrow do
      local guid = dat['guid'][i]
      local tstamp, changeset = sqlConn:rowexec(
	 string.format("select tstamp, GKYL_HG_CHANGESET from RegressionMeta where guid='%s'", guid)
      )
      print(string.format(fmt1, tstamp, changeset, dat['runtime'][i]))
   end
   
end

-- "s" (summary) command
parser:command("summary", "Print summary of all tests")
   :action(summary_action)

-- "q" (query) command
local c_query = parser:command("query", "Query individual run of regression system and print information")
   :action(query_action)
c_query:option("-i --id", "Print information for regression run with this ID", 1)
c_query:option("-t --test", "Print log for specified test number", 0)

-- "h" (history) command
local c_history = parser:command("history", "Query historical data for a specific test")
   :action(history_action)
c_history:option("-r --regression", "Print history for specific test")

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS)
