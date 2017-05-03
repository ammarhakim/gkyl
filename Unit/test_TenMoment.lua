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
   local tenMoment = HyperEquation.TenMoment { }

   assert_equal(10, tenMoment:numEquations(), "No of equations")
   assert_equal(5, tenMoment:numWaves(), "No of waves")

   local rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz = 1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   local q = {rho, rho*u, rho*v, rho*w, pxx+rho*u*u, pxy+rho*u*v, pxz+rho*u*w, pyy+rho*v*v, pyz+rho*v*w, pzz+rho*w*w}
   local flux = Lin.Vec(10)
   tenMoment:flux(1, q, flux)

   local myFlux = {
      rho*u,
      rho*u^2 + pxx,
      rho*u*v + pxy,
      rho*u*w + pxz,
      rho*u^3 + 3*u*pxx,
      rho*u^2*v + 2*u*pxy + v*pxx,
      rho*u^2*w + 2*u*pxz + w*pxx,
      rho*u*v^2 + u*pyy + 2*v*pxy,
      rho*u*v*w + u*pyz + v*pxz + w*pxy,
      rho*u*w^2 + u*pzz + 2*w*pxz
   }

   -- check fluxes
   for i = 1, 10 do
      assert_equal(myFlux[i], flux[i], "Checking flux")
   end

   local rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz = 1.0, 0.5, 0.6, 0.7, 1.0, 0.1, 0.2, 2.0, 0.1, 3.0
   local q = {rho, rho*u, rho*v, rho*w, pxx+rho*u*u, pxy+rho*u*v, pxz+rho*u*w, pyy+rho*v*v, pyz+rho*v*w, pzz+rho*w*w}
   local flux = Lin.Vec(10)
   tenMoment:flux(1, q, flux)

   local myFlux = {
      rho*u,
      rho*u^2 + pxx,
      rho*u*v + pxy,
      rho*u*w + pxz,
      rho*u^3 + 3*u*pxx,
      rho*u^2*v + 2*u*pxy + v*pxx,
      rho*u^2*w + 2*u*pxz + w*pxx,
      rho*u*v^2 + u*pyy + 2*v*pxy,
      rho*u*v*w + u*pyz + v*pxz + w*pxy,
      rho*u*w^2 + u*pzz + 2*w*pxz
   }

   -- check fluxes
   for i = 1, 10 do
      assert_equal(myFlux[i], flux[i], "Checking flux")
   end   
   
end

function test_2()
   
end

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
