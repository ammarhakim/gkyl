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

local function calcq(v)
   local rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz = v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10]
   q = Lin.Vec(10)
   q[1] = rho
   q[2], q[3], q[4] = rho*u, rho*v, rho*w
   q[5] = pxx + rho*u*u
   q[6] = pxy + rho*u*v
   q[7] = pxz + rho*u*w
   q[8] = pyy + rho*v*v
   q[9] = pyz + rho*v*w
   q[10] = pzz + rho*w*w
   return q
end

function test_1()
   local tenMoment = HyperEquation.TenMoment { }

   assert_equal(10, tenMoment:numEquations(), "No of equations")
   assert_equal(5, tenMoment:numWaves(), "No of waves")

   local rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz = 1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   local q = calcq({rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz})
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
   q = calcq({rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz})
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
   local tenMoment = HyperEquation.TenMoment { }

   local ql = calcq({1.0, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0})
   local qr = calcq({1.0, 0.1, 0.1, 0.1, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0})

   local meqn, mwaves = tenMoment:numEquations(), tenMoment:numWaves()
   
   local delta = Lin.Vec(meqn)
   for m = 1, meqn do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(mwaves, meqn)
   local s = Lin.Vec(meqn)
   tenMoment:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(meqn), Lin.Vec(meqn)
   tenMoment:qFluctuations(1, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(meqn)
   for m = 1, meqn do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(meqn), Lin.Vec(meqn), Lin.Vec(meqn)
   tenMoment:flux(1, ql, fluxl)
   tenMoment:flux(1, qr, fluxr)
   for m = 1, meqn do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, meqn do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end   
end

function test_3()
   local tenMoment = HyperEquation.TenMoment { }

   local ql = calcq({1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0})
   local qr = calcq({2.0, 2.0, 0.1, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0})

   local meqn, mwaves = tenMoment:numEquations(), tenMoment:numWaves()
   
   local delta = Lin.Vec(meqn)
   for m = 1, meqn do delta[m] = qr[m]-ql[m] end

   local waves = Lin.Mat(mwaves, meqn)
   local s = Lin.Vec(meqn)
   tenMoment:rp(1, delta, ql, qr, waves, s)
   local amdq, apdq = Lin.Vec(meqn), Lin.Vec(meqn)
   tenMoment:qFluctuations(1, ql, qr, waves, s, amdq, apdq)

   local sumFluct = Lin.Vec(meqn)
   for m = 1, meqn do sumFluct[m] = amdq[m]+apdq[m] end
   local fluxr, fluxl, df = Lin.Vec(meqn), Lin.Vec(meqn), Lin.Vec(meqn)
   tenMoment:flux(1, ql, fluxl)
   tenMoment:flux(1, qr, fluxr)
   for m = 1, meqn do df[m] = fluxr[m]-fluxl[m] end

   for m = 1, meqn do
      assert_equal(df[m], sumFluct[m], "Checking jump in flux is sum of fluctuations")
   end   
end

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
