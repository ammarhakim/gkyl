-- Gkyl ------------------------------------------------------------------------
--
-- Test for Advection equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Advection = require "Eq.Advection"
local Basis = require "Basis"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local advection = Advection { ndim = 2, velocity = {1.0, 2.0} }

   assert_equal(1, advection:numEquations(), "Testing number of equations")
   assert_equal(1, advection:numWaves(), "Testing number of waves")

   -- dir 1
   local ql, qr = {1.0}, {2.5}
   local delta, s, waves = Lin.Vec(1), Lin.Vec(1), Lin.Mat(1, 1)
   delta[1] = qr[1]-ql[1]
   advection:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(1), Lin.Vec(1)
   advection:qFluctuations(1, ql, q2, waves, s, amdq, apdq)

   local fluxl, fluxr, df = Lin.Vec(1), Lin.Vec(1), Lin.Vec(1)
   advection:flux(1, ql, fluxl)
   advection:flux(1, qr, fluxr)
   df[1] = fluxr[1]-fluxl[1]

   assert_equal(df[1], apdq[1]+amdq[1], "Checking jump in flux is sum of fluctuations")

   -- dir 2
   advection:rp(2, delta, ql, qr, waves, s)
   advection:qFluctuations(1, ql, q2, waves, s, amdq, apdq)

   advection:flux(2, ql, fluxl)
   advection:flux(2, qr, fluxr)
   df[1] = fluxr[1]-fluxl[1]

   assert_equal(df[1], apdq[1]+amdq[1], "Checking jump in flux is sum of fluctuations")   
end

function test_2()
   local advection = Advection { ndim = 2, velocity = {-1.0, -2.0} }

   -- dir 1
   local ql, qr = {1.0}, {2.5}
   local delta, s, waves = Lin.Vec(1), Lin.Vec(1), Lin.Mat(1, 1)
   delta[1] = qr[1]-ql[1]
   advection:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(1), Lin.Vec(1)
   advection:qFluctuations(1, ql, q2, waves, s, amdq, apdq)

   local fluxl, fluxr, df = Lin.Vec(1), Lin.Vec(1), Lin.Vec(1)
   advection:flux(1, ql, fluxl)
   advection:flux(1, qr, fluxr)
   df[1] = fluxr[1]-fluxl[1]

   assert_equal(df[1], apdq[1]+amdq[1], "Checking jump in flux is sum of fluctuations")

   -- dir 2
   advection:rp(2, delta, ql, qr, waves, s)
   advection:qFluctuations(1, ql, q2, waves, s, amdq, apdq)

   advection:flux(2, ql, fluxl)
   advection:flux(2, qr, fluxr)
   df[1] = fluxr[1]-fluxl[1]

   assert_equal(df[1], apdq[1]+amdq[1], "Checking jump in flux is sum of fluctuations")   
end

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
