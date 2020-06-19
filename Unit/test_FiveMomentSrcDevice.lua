-- Gkyl ------------------------------------------------------------------------
--
-- Test for the gpu implementation of the time-centered five-moment source
-- updater against the cpu implementation.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local cuda
local cuAlloc
if GKYL_HAVE_CUDA then
  cuda = require "Cuda.RunTime"
  cuAlloc = require "Cuda.Alloc"
end

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Time = require "Lib.Time"
local xsys = require "xsys"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

local nx = 128*64 -- number of configuration space dimensions in x
local nloop = NLOOP or 3 -- number of WavePropagation calls to loop over
local numThreads = NTHREADS or 128 -- number of threads to use in WavePropagation kernel configuration
local checkResult = xsys.pickBool(CHECK, true) and false

local grid = Grid.RectCart {
   lower = {0.0},
   upper = {1.0},
   cells = {nx},
}

local fluid1 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
}
local fluid2 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
}
local emf = DataStruct.Field {
   onGrid = grid,
   numComponents = 8,
}

local d_fluid1 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   createDeviceCopy = true,
}
local d_fluid2 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   createDeviceCopy = true,
}
local d_emf = DataStruct.Field {
   onGrid = grid,
   numComponents = 8,
   createDeviceCopy = true,
}

local epsilon0 = 1.0
local gasGamma = 5. / 3.

local computeMaxFreq = function(charge, mass, rho, Bmag, epsilon0, show)
   local wp_tot = 0
   local wc_max = 0
   for species,_ in ipairs(charge)  do
      local wp2 = rho[species] * (charge[species] / mass[species])^2 / epsilon0
      wp_tot = wp_tot + wp2
      wc_max = math.max(wc_max, math.abs(charge[species] * Bmag / mass[species]))
   end
   wp_tot = math.sqrt(wp_tot)

   if show then
      for species = 1,#charge do
         print(string.format("species [%d] charge=%g, mass=%g, rho=%g",
            species, charge[species], mass[species], rho[species]))
      end
      print('Bmag=%g, epsilon0=%g')

      print('wp_tot', wp_tot)
      print('2*pi/wp_tot', 2*math.pi/wp_tot)
      print('wc_max', wc_max)
      print('2*pi/wc_max', 2*math.pi/wc_max)
   end

   return math.max(wp_tot, wc_max)
end

function init(
   fluid1, fluid2, emf,
   rho1, rhovx1, rhovy1, rhovz1, p1,
   rho2, rhovx2, rhovy2, rhovz2, p2,
   Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
   )

   local f1Idxr = fluid1:genIndexer()   
   local f2Idxr = fluid2:genIndexer()   
   local emIdxr = emf:genIndexer() 
   for idx in fluid1:localRangeIter() do
      local f1 = fluid1:get(f1Idxr(idx))
      local f2 = fluid2:get(f2Idxr(idx))
      local em = emf:get(emIdxr(idx))

      local e1 = 0.5 * (rhovx1^2 + rhovy1^2 + rhovz1^2) / rho1 + p1 / (gasGamma - 1)
      f1[1] = rho1
      f1[2] = rhovx1
      f1[3] = rhovy1
      f1[4] = rhovz1
      f1[5] = e1

      local e2 = 0.5 * (rhovx2^2 + rhovy2^2 + rhovz2^2) / rho2 + p2 / (gasGamma - 1)
      f2[1] = rho2
      f2[2] = rhovx2
      f2[3] = rhovy2
      f2[4] = rhovz2
      f2[5] = e2

      em[1] = Ex
      em[2] = Ey
      em[3] = Ez
      em[4] = Bx
      em[5] = By
      em[6] = Bz
      em[7] = phiE
      em[8] = phiB
   end
   fluid1:copyHostToDevice()
   fluid2:copyHostToDevice()
   emf:copyHostToDevice()

end

function doTest1(scheme, dtFrac)

   local charge1, mass1, rho1, rhovx1, rhovy1, rhovz1, p1 = -1, 1, 4, 1, 2, 3, 1
   local charge2, mass2, rho2, rhovx2, rhovy2, rhovz2, p2 = 1, 4, 16, 6, 8, 10, 1.5
   local Ex, Ey, Ez = 1.2, 1.3, 1.4
   local Bx, By, Bz = 1, 2, 3
   local phiE, phiB = 3.33, 4.56
   local Bmag = math.sqrt(Bx^2 + By^2 + Bz^2)

   local charge = {charge1, charge2}
   local mass ={mass1, mass2}
   local rho = {rho1, rho2}
   local wmax = computeMaxFreq(charge, mass, rho, Bmag, epsilon0)
   local dt = dtFrac/wmax

   init(
      fluid1, fluid2, emf,
      rho1, rhovx1, rhovy1, rhovz1, p1,
      rho2, rhovx2, rhovy2, rhovz2, p2,
      Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
   )

   init(
      d_fluid1, d_fluid2, d_emf,
      rho1, rhovx1, rhovy1, rhovz1, p1,
      rho2, rhovx2, rhovy2, rhovz2, p2,
      Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
   )
   d_fluid1:copyHostToDevice()
   d_fluid2:copyHostToDevice()
   d_emf:copyHostToDevice()

   local srcUpdater = Updater.FiveMomentSrc {
      onGrid = grid,
      numFluids = #charge,
      charge = charge,
      mass = mass,
      epsilon0 = epsilon0,
      elcErrorSpeedFactor = 1,
      mgnErrorSpeedFactor = 1,
      scheme = scheme,
   }
   srcUpdater:setDtAndCflRate(dt, nil)

   tmStart = Time.clock()
   for i = 1, nloop do
      srcUpdater:_advanceDevice(0.0, {}, {d_fluid1, d_fluid2, d_emf})
   end
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)

   d_fluid1:copyDeviceToHost()
   d_fluid2:copyDeviceToHost()
   d_emf:copyDeviceToHost()

   tmStart = Time.clock()
   for i = 1, nloop do
      srcUpdater:advance(0.0, {}, {fluid1, fluid2, emf})
   end
   local totalCpuTime = (Time.clock()-tmStart)

   if checkResult then
      local f1Idxr = fluid1:genIndexer()   
      local f2Idxr = fluid2:genIndexer()   
      local emIdxr = emf:genIndexer() 
      local d_f1Idxr = d_fluid1:genIndexer()   
      local d_f2Idxr = d_fluid2:genIndexer()   
      local d_emIdxr = d_emf:genIndexer() 

      for idx in fluid1:localRangeIter() do
         local f1 = fluid1:get(f1Idxr(idx))
         local d_f1 = d_fluid1:get(d_f1Idxr(idx))
         for comp = 1,5 do
            assert_close(f1[comp], d_f1[comp], 1e-10, string.format(
            "fluid1 index %d component %d is incorrect", f1Idxr(idx), comp))
         end

         local f2 = fluid1:get(f2Idxr(idx))
         local d_f2 = d_fluid1:get(d_f2Idxr(idx))
         for comp = 1,5 do
            assert_close(f2[comp], d_f2[comp], 1e-10, string.format(
            "fluid2 index %d component %d is incorrect", f2Idxr(idx), comp))
         end

         local em = fluid1:get(emIdxr(idx))
         local d_em = d_fluid1:get(d_emIdxr(idx))
         for comp = 1,8 do
            assert_close(em[comp], d_em[comp], 1e-10, string.format(
            "emf index %d component %d is incorrect", emIdxr(idx), comp))
         end
      end
   end

   print("cpu scheme", scheme)
   print(string.format("Total CPU time for %d FiveMomentSrc calls = %f s   (average = %f s)", nx, totalCpuTime, totalCpuTime/nx))
   print(string.format("Total GPU time for %d FiveMomentSrc calls = %f s   (average = %f s)", nx, totalGpuTime, totalGpuTime/nx))
   print(string.format("GPU speed-up = %fx!", totalCpuTime/totalGpuTime))

end

doTest1("time-centered", 1)
doTest1("direct", 1)

-- TODO Do automatic check.
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
