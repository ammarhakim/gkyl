-- Gkyl ------------------------------------------------------------------------
--
-- Debugging script to dump tables created by runregression code
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local sql = require "sqlite3"

local f = loadfile("runregression.config.lua")
if not f then
   log("Regression tests not configured! Run config command first\n")
   os.exit(1)
end
configVals = f()

-- open connection
sqlConn = sql.open(string.format("%s/regressiondb", configVals['results_dir']))

sqlConn "select * from RegressionMeta"
