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

-- function to handle summary ("s") command
local function summary_action(args, name)
   local dbMeta = read_metatable()
   local nrow = #dbMeta
      
   print("---------------------------------------------------")
   print("ID: Time-Stamp Changeset Build-Date Total Pass Fail")
   print("---------------------------------------------------")
   for i, d in pairs(dbMeta) do
      --print(d.tstamp, d.changeset, d.builddate, d.ntotal, d.npass, d.nfail)
      print(string.format("%d: %s %s %s %s %s %s",
        nrow-i+1, d.tstamp, d.changeset, d.builddate, tostring(d.ntotal), tostring(d.npass), tostring(d.nfail)))
   end
end

-- function to handle query ("q") command
local function query_action(args, name)
   local dbMeta = read_metatable()
   local id = #dbMeta-args.id+1 -- meta-list is printed in reverse order
   local guid = dbMeta[id].guid -- this is tha GUID for fetching information

   local t, nrow = sqlConn:exec(
      string.format("select * from RegressionData where guid=='%s'", guid))

   local function status2str(status)
      local s = tonumber(status)
      if s == -2 then return "create" end
      if s == -1 then return "skip  " end      
      if s == 0 then return  "fail  " end
      if s == 1 then return  "pass  " end
   end

   
   print("------------------------")
   print("ID: Name Status Run-Time")
   print("------------------------")   
   for i = 1, nrow do
      print(string.format(
	       "%d: %s %s %g sec", i, t['name'][i], status2str(t['status'][i]), t['runtime'][i]))
   end
end

-- "s" (summary) command
parser:command("s", "Print summary of all tests")
   :action(summary_action)

-- "q" (query) command
local c_query = parser:command("q", "Query individual run of regression system and print information")
   :action(query_action)
c_query:option("-i --id", "Print information for regression run with this ID")

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS)
