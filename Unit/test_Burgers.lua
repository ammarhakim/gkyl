-- Gkyl ------------------------------------------------------------------------
--
-- Test for Burgers equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local HyperEquation = require "Eq"
local Basis = require "Basis"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local burgers = HyperEquation.Burgers {}

   assert_equal(1, burgers:numEquations(), "Testing number of equations")
   assert_equal(1, burgers:numWaves(), "Testing number of waves")

   local ql, qr = {1.0}, {2.5}
   local delta, s, waves = Lin.Vec(1), Lin.Vec(1), Lin.Mat(1, 1)
   delta[1] = qr[1]-ql[1]
   burgers:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(1), Lin.Vec(1)
   burgers:qFluctuations(1, ql, q2, waves, s, amdq, apdq)   

   local fluxl, fluxr, df = Lin.Vec(1), Lin.Vec(1), Lin.Vec(1)
   burgers:flux(1, ql, fluxl)
   burgers:flux(1, qr, fluxr)
   df[1] = fluxr[1]-fluxl[1]

   assert_equal(df[1], apdq[1]+amdq[1], "Checking jump in flux is sum of fluctuations")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
