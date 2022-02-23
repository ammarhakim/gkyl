-- Gkyl ------------------------------------------------------------------------
--
-- Test for Vlasov equation object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Basis = require "Basis"
local ffi  = require "ffi"
local Unit = require "Unit"
local Vlasov = require "Eq.Vlasov"
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local confBasis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local phaseBasis = Basis.CartModalSerendipity { ndim = 3, polyOrder = 2 }
   local confExtRange = Range.Range({0}, {10})
   local confRange = confExtRange:subRange({1}, {9})
   
   local vlasovEq = Vlasov {
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      confRange = confRange,
      charge = -1.0,
      mass = 1.0,
   }

   assert_equal(1, vlasovEq:numEquations(), "Checking number of equations")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
