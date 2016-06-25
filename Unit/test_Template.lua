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

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
