-- Gkyl ------------------------------------------------------------------------
--
-- Test for Perfectly hyperbolic Maxwell equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

if GKYL_HAVE_CUDA == false then
   print("**** Can't run kernel tests without CUDA enabled GPUs!")
   return 0
end

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

   local ql = Lin.Vec(8)
   local qr = Lin.Vec(8)
   for m = 1, 8 do
      ql[m] = math.random()
      qr[m] = math.random()
   end
   -- local _ql = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}
   -- local _qr = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}
   -- for m = 1, 8 do
   --    ql[m] = _ql[m]
   --    qr[m] = _qr[m]
   -- end

   local delta = Lin.Vec(8)
   for m = 1, 8 do delta[m] = qr[m]-ql[m] end

   -- C vs. Lua: rp and qFluctuations
   local waves = Lin.Mat(6, 8)
   local s = Lin.Vec(8)
   maxwell:rpCImpl(0, delta, ql._p, qr._p, waves._p, s._p)

   local waves0 = Lin.Mat(6, 8)
   local s0 = Lin.Vec(8)
   maxwell:rp(1, delta, ql, qr, waves0, s0)

   for w = 1, 6 do
      for m = 1, 8 do
         assert_equal(waves0[w][m], waves[w][m],
            string.format("Checking waves [%d][%d]", w, m))
      end
   end

   for w = 1, 6 do
      assert_equal(s0[w], s[w], string.format("Checking speeds [%d]", w))
   end

   local amdq, apdq = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuationsCImpl(0, ql._p, qr._p, waves._p, s._p, amdq._p, apdq._p)
   local amdq0, apdq0 = Lin.Vec(8), Lin.Vec(8)
   maxwell:qFluctuations(1, ql, qr, waves0, s0, amdq0, apdq0)

   for m = 1, 8 do
      assert_equal(amdq0[m], amdq[m], string.format("Checking fluctuation amdq [%d]", m))
   end

   for m = 1, 8 do
      assert_equal(apdq0[m], apdq[m], string.format("Checking fluctuation apdq [%d]", m))
   end

   local sumFluct = Lin.Vec(8)
   for m = 1, 8 do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(8), Lin.Vec(8), Lin.Vec(8)
   maxwell:fluxCImpl(0, ql._p, fluxl._p)
   maxwell:fluxCImpl(0, qr._p, fluxr._p)
   for m = 1, 8 do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, 8 do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end   
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
