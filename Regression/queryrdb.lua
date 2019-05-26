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

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS)
