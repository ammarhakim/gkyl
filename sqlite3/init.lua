-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the sqlite3 library module
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

--[[

 See: http://scilua.org/ljsqlite3.html for API documentation. Except
 one needs to load this as

 local sql = require "sqlite3"
--]]

return require("sqlite3.sqlite3")
