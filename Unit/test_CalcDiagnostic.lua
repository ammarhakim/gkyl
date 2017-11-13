-- Gkyl ------------------------------------------------------------------------
--
-- Test for diagnostic updater
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CalcDiagnostic = require "Updater.CalcDiagnostic"
local ffi  = require "ffi"

local assert_equal = Unit.assert_equal
local stats = Unit.stats
