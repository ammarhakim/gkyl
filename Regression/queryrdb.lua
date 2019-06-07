-- Gkyl ------------------------------------------------------------------------
--
-- Query code to get information from DB created by runregression code
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local sql = require "sqlite3"

-- SQLite connection handle
local sqlConn = nil
-- status code -> string
local statusToString = { [-2] = "create", [-1] = "skip", [0] = "fail", [1] = "pass" }

-- configure query system
local function configure(args)
   if args.db then
      sqlConn = sql.open(args.db)
   else
      -- load configuration file created by regression testing system
      local f = loadfile("runregression.config.lua")
      if not f then
	 log("Regression tests not configured! Run config command first\n")
	 os.exit(1)
      end
      local configVals = f()
      sqlConn = sql.open(string.format("%s/regressiondb", configVals['results_dir']))
   end
end

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
   local maxNm = 0
   dbData = {}
   for i = 1, nrow do
      dbData[i] = {
	 name = t['name'][i],
	 status = statusToString[tonumber(t['status'][i])],
	 runtime = t['runtime'][i],
	 runlog = t['runlog'][i],
      }
      maxNm = math.max(maxNm, string.len(t['name'][i]))
   end
   return dbData, maxNm
end

-- function to handle summary command
local function summary_action(args, name)
   configure(args)
   
   local dbMeta = read_metatable()
   local nrow = #dbMeta

   local fmt = "%-4s: %-20s %-30s %-5s %-5s %-5s"
   print(string.format(fmt, "ID", "Time-Stamp", "Changeset", "Total", "Pass", "Fail"))
   for i,d in pairs(dbMeta) do
      print(
	 string.format(
	    fmt,
	    nrow-i+1, d.tstamp, d.changeset, tonumber(d.ntotal), tonumber(d.npass), tonumber(d.nfail)
	 )
      )
   end
end

-- function to handle query command
local function query_action(args, name)
   configure(args)
   
   local dbMeta = read_metatable()
   local idx = #dbMeta-args.id+1 -- meta-list is printed in reverse order


   if args.net_time then
      -- print total time to run all tests specified ID
      local t, nrow = sqlConn:exec(
	 string.format("select runtime from RegressionData where guid=='%s'", 
		       dbMeta[idx].guid))
      local nettm = 0.0
      for i = 1, nrow do
	 nettm = nettm+t['runtime'][i]
      end
      print(string.format("%.4g", nettm))
   else
      local dbData, maxNm = read_tests_with_id(dbMeta[idx].guid)
      -- print info about all tests with specified ID
      local function filt(status)
	 if args.fail_only then
	    if status == "fail" then return true else return false end
	 end
	 return true
      end
      
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
	    if filt(d.status) then
	       print(string.format(fmt1, i, d.name, d.status, d.runtime))
	    end
	 end
      end
   end
end

-- function to handle delete command
local function delete_action(args, name)
   configure(args)
   if args.id < 1 then
      print("No deletions done. Specify run ID as returned by summary command")
   end
end

-- function to handle history command
local function history_action(args, name)
   configure(args)
   
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

   local fmt = "%-20s %-30s %-7s %-4s"
   print(string.format(fmt, "Time-Stamp", "Changeset", "Status", "Run-Time"))
   local fmt1 = "%-20s %-30s %-7s %.4g"
   for i = 1,nrow do
      local guid = dat['guid'][i]
      local tstamp, changeset = sqlConn:rowexec(
	 string.format("select tstamp, GKYL_HG_CHANGESET from RegressionMeta where guid='%s'", guid)
      )
      local stat = statusToString[tonumber(dat['status'][i])]
      print(string.format(fmt1, tstamp, changeset, stat, dat['runtime'][i]))
   end
   
end

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("queryrdb")
   :require_command(true)
   :description [[
Query database created by runregression system to extract information
about various tests.
]]
parser:option("--db", "Alternate DB to use")

-- summary command
parser:command("summary", "Print summary of all tests")
   :action(summary_action)

-- query command
local c_query = parser:command("query", "Query individual run of regression system and print information")
   :action(query_action)
c_query:option("-i --id", "Print information for regression run with this ID", 1)
c_query:flag("-f --fail-only", "Print only failed tests", false)
c_query:option("-t --test", "Print log for specified test number", 0)
c_query:flag("--net-time", "Print total time it too to run all tests", false)

-- delete command
local c_delete = parser:command("delete", "Delete run data from regression system")
   :action(delete_action)
c_delete:option("-i --id", "Regression run with this ID will be deleted", 0)

-- history command
local c_history = parser:command("history", "Query historical data for a specific test")
   :action(history_action)
c_history:option("-r --regression", "Print history for specific test")

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS)
