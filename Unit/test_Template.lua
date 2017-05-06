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

local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   --local w1, w2, w3 = waves[1], waves[2], waves[3]
   local s1m, s2m, s3m = math.min(0, s[1]), math.min(0, s[2]), math.min(0, s[3])
   local s1p, s2p, s3p = math.max(0, s[1]), math.max(0, s[2]), math.max(0, s[3])

|for i = 1, 5 do
   amdq[${i}] = s1m*waves[1][${i}] + s2m*waves[2][${i}] + s3m*waves[3][${i}]
   apdq[${i}] = s1p*waves[1][${i}] + s2p*waves[2][${i}] + s3p*waves[3][${i}]

|end
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = loadstring( qFluctuationsTempl {} )()

local intersectTempl = xsys.template([[
return function (s, rgn)
| for d = 1, NDIM do
	 if math.min(s:upper(${d}), rgn:upper(${d})) < math.max(s:lower(${d}), rgn:lower(${d})) then
	    return true
	 end
| end
  return false
end
]])
-- instantiate it
local intersect = intersectTempl {NDIM=3}

local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   local w1, w2, w3, w4, w5 = waves[1], waves[2], waves[3], waves[4], waves[5]
   local s1m, s2m, s3m, s4m, s5m = math.min(0, s[1]), math.min(0, s[2]), math.min(0, s[3]), math.min(0, s[4]), math.min(0, s[5])
   local s1p, s2p, s3p, s4p, s5p = math.max(0, s[1]), math.max(0, s[2]), math.max(0, s[3]), math.max(0, s[4]), math.max(0, s[5])

|for i = 1, 10 do
   amdq[${i}] = s1m*w1[${i}] + s2m*w2[${i}] + s3m*w3[${i}] + s4m*w4[${i}] + s5m*w5[${i}]
   apdq[${i}] = s1p*w1[${i}] + s2p*w2[${i}] + s3p*w3[${i}] + s4p*w4[${i}] + s5p*w5[${i}]

|end
end
]])
local qFluctuations = qFluctuationsTempl {}

local waveDotProdTempl = xsys.template([[
return function (meqn, waves, waves1, mw)
  local mw1 = mw-1
  return
|for i = 0, MEQN-2 do
  waves1[${MEQN}*mw1+${i}]*waves[${MEQN}*mw1+${i}]+
|end
  waves1[${MEQN}*mw1+${MEQN-1}]*waves[${MEQN}*mw1+${MEQN-1}]
end
]])
local waveDotProd = waveDotProdTempl { MEQN=5 }


-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   local w1, w2, w3, w4, w5 = waves[1], waves[2], waves[3], waves[4], waves[5]
   local s1m, s2m, s3m, s4m, s5m = math.min(0, s[1]), math.min(0, s[2]), math.min(0, s[3]), math.min(0, s[4]), math.min(0, s[5])
   local s1p, s2p, s3p, s4p, s5p = math.max(0, s[1]), math.max(0, s[2]), math.max(0, s[3]), math.max(0, s[4]), math.max(0, s[5])

|for i = 1, 10 do
   amdq[${i}] = s1m*w1[${i}] + s2m*w2[${i}] + s3m*w3[${i}] + s4m*w4[${i}] + s5m*w5[${i}]
   apdq[${i}] = s1p*w1[${i}] + s2p*w2[${i}] + s3p*w3[${i}] + s4p*w4[${i}] + s5p*w5[${i}]

|end
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = qFluctuationsTempl {}
print(qFluctuations)

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
