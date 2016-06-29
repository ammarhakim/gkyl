-- Gkyl ------------------------------------------------------------------------
--
-- Test for templating library
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local Unit = require "Unit"
local HyperEquation = require "Eq"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

local calcDeltaTempl = xsys.template([[
return function (ql, qr, delta)
|for i = 1, MEQN do
  delta[${i}] = qr[${i}] - ql[${i}]
|end
end
]])
function test_1()
   local calcDelta = loadstring( calcDeltaTempl {MEQN = 5})()
   local ql, qr, delta = Lin.Vec(5), Lin.Vec(5), Lin.Vec(5)

   ql[1], ql[2], ql[3], ql[4], ql[5] = 0, 0, 0, 0, 0
   qr[1], qr[2], qr[3], qr[4], qr[5] = 1, 2, 3, 4, 5
   calcDelta(ql, qr, delta)

   assert_equal(1, delta[1], "Testing calcDelta")
   assert_equal(2, delta[2], "Testing calcDelta")
   assert_equal(3, delta[3], "Testing calcDelta")
   assert_equal(4, delta[4], "Testing calcDelta")
   assert_equal(5, delta[5], "Testing calcDelta")
end

local calcCflaTempl = xsys.template([[
return function (cfla, dtdx, s)
  local c = cfla
|for i = 1, MWAVE do
  c = math.max(c, dtdx*math.abs(s[${i}]))
|end
  return c
end
]])
function test_2()
   local calcCflaTempl = loadstring( calcCflaTempl {MWAVE = 3})()
   local s = Lin.Vec(3)

   s[1], s[2], s[3] = -2, 0.1, 3
   local cfla = calcCflaTempl(0, 1.0, s)
   assert_equal(3, cfla, "Testing calcCflaTempl")

   s[1], s[2], s[3] = -3.5, 0.1, 3
   local cfla = calcCflaTempl(0, 1.0, s)
   assert_equal(3.5, cfla, "Testing calcCflaTempl")

   s[1], s[2], s[3] = -2, 0.1, 3
   local cfla = calcCflaTempl(10, 1.0, s)
   assert_equal(10, cfla, "Testing calcCflaTempl")
end

local calcFirstOrderGudTempl = xsys.template([[
return function (dtdx, ql, qr, amdq, apdq)
|for i = 1, MEQN do
  qr[${i}] = qr[${i}] - dtdx*apdq[${i}]
|end
|for i = 1, MEQN do
  ql[${i}] = ql[${i}] - dtdx*amdq[${i}]
|end
end
]])

local waveDotProdTempl = xsys.template([[
return function (meqn, waves, waves1, mw)
  local mw1 = mw-1
  return
|for i = 0, MEQN-2 do
  waves[meqn*mw1+${i}]*waves[meqn*mw1+${i}]+
|end
  waves1[meqn*mw1+${MEQN-1}]*waves[meqn*mw1+${MEQN-1}]
end
]])

function test_3()
   local mwave, meqn = 3, 5
   local waveDotProd = loadstring( waveDotProdTempl {MEQN=meqn} )()
   local waves = ffi.new("double[?]", mwave*meqn)

   for i = 1, mwave*meqn do waves[i-1] = i end

   local dotr = {1^2+2^2+3^2+4^2+5^2,
		 6^2+7^2+8^2+9^2+10^2,
		 11^2+12^2+13^2+14^2+15^2}
   for mw = 1, 3 do
      assert_equal(dotr[mw], waveDotProd(meqn, waves, waves, mw), string.format("Testing dot product of wave %d", mw))
   end
end

local rescaleWaveTempl = xsys.template([[
return function (scale, wave)
|for i = 0, MEQN-1 do
  wave[${i}] = scale*wave[${i}]
|end
end
]])

local secondOrderTempl = xsys.template([[
return function (dtdx, s, wave, fs)
  local sfact = 0.5*math.abs(s)*(1-math.abs(s)*dtdx)
|for i = 0, MEQN-1 do
  fs[${i}] = fs[${i}] + sfact*wave[${i}]
|end
end
]])

local rotateWavesToGlobalTempl = xsys.template([[
return function (dir, equation, wavesLocal, waves)
|for i = 1, MWAVE do
  local wIn, wOut = wavesLocal[${i}], waves[${i}]
  equation:rotateToGlobalAligned(dir, wIn, wOut)

|end
end
]])

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end

print(rotateWavesToGlobalTempl {MWAVE = 3})
