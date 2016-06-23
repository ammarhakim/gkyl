-- Gkyl ------------------------------------------------------------------------
--
-- Test for Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local HyperEquation = require "Eq"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local euler = HyperEquation.Euler { gasGamma = 1.4 }

   assert_equal(5, euler:numEquations(), "No of equations")
   assert_equal(3, euler:numWaves(), "No of wave")

   local rho, u, v, w, pr = 1.0, 1.0, 2.0, 3.0, 0.1
   local q = Lin.Vec(5)
   q[1] = rho
   q[2], q[3], q[4] = rho*u, rho*v, rho*w
   q[5] = pr/(euler:gasGamma()-1) + 0.5*rho*(u*u+v*v+w*w)

   local prOut = euler:pressure(q)
   assert_equal(pr, prOut, "Checking pressure")

   local flux = Lin.Vec(5)
   euler:flux(q, flux)
   assert_equal(flux[1], rho*u, "Checking flux")
   assert_equal(flux[2], rho*u*u+pr, "Checking flux")
   assert_equal(flux[3], rho*u*v, "Checking flux")
   assert_equal(flux[4], rho*u*w, "Checking flux")
   assert_equal(flux[5], (q[5]+pr)*u, "Checking flux")

   rho , u, v, w, pr = -1.0, 1.0, 0.0, 0.0, 0.1
   q[1] = rho
   q[2], q[3], q[4] = rho*u, rho*v, rho*w
   q[5] = pr/(euler:gasGamma()-1) + 0.5*rho*(u*u+v*v+w*w)
   assert_equal(false, euler:isPositive(q), "Checking positivity")

   rho , u, v, w, pr = 1.0, 1.0, 0.0, 0.0, -0.1
   q[1] = rho
   q[2], q[3], q[4] = rho*u, rho*v, rho*w
   q[5] = pr/(euler:gasGamma()-1) + 0.5*rho*(u*u+v*v+w*w)
   assert_equal(false, euler:isPositive(q), "Checking positivity")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
