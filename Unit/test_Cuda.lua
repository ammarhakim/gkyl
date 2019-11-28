-- Gkyl ------------------------------------------------------------------------
--
-- Test for CUDA runtime API wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't do anything if we were not built with CUDA
if GKYL_HAVE_CUDA == false then
    return 0
end

local Unit = require "Unit"
local cuda = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
