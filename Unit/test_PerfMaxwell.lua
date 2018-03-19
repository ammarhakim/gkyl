-- Gkyl ------------------------------------------------------------------------
--
-- Test for Perfectly hyperbolic Maxwell equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local PerfMaxwell = require "Eq.PerfMaxwell"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local maxwell = PerfMaxwell {
      lightSpeed = 5.0,
      elcErrorSpeedFactor = 1.0,
      mgnErrorSpeedFactor = 3.0,
   }

   local ql = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 4.0}
   local qr = {0.1, 0.1, 0.0, 0.0, 2.0, 3.0, 3.0, 3.0}

   local delta = Lin.Vec(8)
   for m = 1, 8 do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(6, 8)
   local s = Lin.Vec(8)
   maxwell:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuations(1, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(8)
   for m = 1, 8 do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(8), Lin.Vec(8), Lin.Vec(8)
   maxwell:flux(1, ql, fluxl)
   maxwell:flux(1, qr, fluxr)
   for m = 1, 8 do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, 8 do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end   
end

function test_2()
   local maxwell = PerfMaxwell {
      lightSpeed = 5.0,
      elcErrorSpeedFactor = 1.0,
      mgnErrorSpeedFactor = 3.0,
   }

   local ql = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}
   local qr = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}

   local delta = Lin.Vec(8)
   for m = 1, 8 do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(6, 8)
   local s = Lin.Vec(8)
   maxwell:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuations(1, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(8)
   for m = 1, 8 do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(8), Lin.Vec(8), Lin.Vec(8)
   maxwell:flux(1, ql, fluxl)
   maxwell:flux(1, qr, fluxr)
   for m = 1, 8 do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, 8 do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end   
end

function test_3()
   local maxwell = PerfMaxwell {
      lightSpeed = 5.0,
      elcErrorSpeedFactor = 1.0,
      mgnErrorSpeedFactor = 3.0,
   }

   local ql = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}
   local qr = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}

   local delta = Lin.Vec(8)
   for m = 1, 8 do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(6, 8)
   local s = Lin.Vec(8)
   maxwell:rp(2, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuations(2, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(8)
   for m = 1, 8 do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(8), Lin.Vec(8), Lin.Vec(8)
   maxwell:flux(2, ql, fluxl)
   maxwell:flux(2, qr, fluxr)
   for m = 1, 8 do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, 8 do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end
end

function test_4()
   local maxwell = PerfMaxwell {
      lightSpeed = 5.0,
      elcErrorSpeedFactor = 1.0,
      mgnErrorSpeedFactor = 3.0,
   }

   local ql = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}
   local qr = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}

   local delta = Lin.Vec(8)
   for m = 1, 8 do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(6, 8)
   local s = Lin.Vec(8)
   maxwell:rp(3, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuations(3, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(8)
   for m = 1, 8 do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(8), Lin.Vec(8), Lin.Vec(8)
   maxwell:flux(3, ql, fluxl)
   maxwell:flux(3, qr, fluxr)
   for m = 1, 8 do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, 8 do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end
end

test_1()
test_2()
test_3()
test_4()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
