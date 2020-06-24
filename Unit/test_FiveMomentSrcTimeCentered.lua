-- Gkyl ------------------------------------------------------------------------
--
-- Test for the different cpu implementations of the time-centered
-- five-moment source updater.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Time = require "Lib.Time"

local assert_close = Unit.assert_close
local stats = Unit.stats

local nloop = NLOOP or 1 -- number of WavePropagation calls to loop over

local grid = Grid.RectCart {
   lower = {0.0},
   upper = {1.0},
   cells = {1},
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
}
local d_fluid2 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
}
local d_emf = DataStruct.Field {
   onGrid = grid,
   numComponents = 8,
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

function doTest1(scheme, d_scheme, dtFrac)

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

   local d_srcUpdater = Updater.FiveMomentSrc {
      onGrid = grid,
      numFluids = #charge,
      charge = charge,
      mass = mass,
      epsilon0 = epsilon0,
      elcErrorSpeedFactor = 1,
      mgnErrorSpeedFactor = 1,
      scheme = d_scheme,
   }
   d_srcUpdater:setDtAndCflRate(dt, nil)

   tmStart = Time.clock()
   for i = 1, nloop do
      d_srcUpdater:_advance(0.0, {}, {d_fluid1, d_fluid2, d_emf})
   end
   local totalGpuTime = (Time.clock()-tmStart)

   tmStart = Time.clock()
   for i = 1, nloop do
      srcUpdater:advance(0.0, {}, {fluid1, fluid2, emf})
   end
   local totalCpuTime = (Time.clock()-tmStart)

   local good = 0
   local bad = 0
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
         if math.abs(f1[comp]-d_f1[comp])>1e-10 then
            bad = bad + 1
         else
            good = good + 1
         end
         assert_close(f1[comp], d_f1[comp], 1e-10, string.format(
         "fluid1 index %d component %d is incorrect", f1Idxr(idx), comp))
      end

      local f2 = fluid1:get(f2Idxr(idx))
      local d_f2 = d_fluid1:get(d_f2Idxr(idx))
      for comp = 1,5 do
         if math.abs(f2[comp]-d_f2[comp])>1e-10 then
            bad = bad + 1
         else
            good = good + 1
         end
         assert_close(f2[comp], d_f2[comp], 1e-10, string.format(
         "fluid2 index %d component %d is incorrect", f2Idxr(idx), comp))
      end

      local em = fluid1:get(emIdxr(idx))
      local d_em = d_fluid1:get(d_emIdxr(idx))
      for comp = 1,8 do
         if math.abs(em[comp]-d_em[comp])>1e-10 then
            bad = bad + 1
         else
            good = good + 1
         end
         assert_close(em[comp], d_em[comp], 1e-10, string.format(
         "emf index %d component %d is incorrect", emIdxr(idx), comp))
      end
   end
   print("good values", good, "bad values", bad)
   print("scheme:", srcUpdater.scheme, "vs", d_srcUpdater.scheme)

end

-- doTest1("time-centered", 1)
doTest1("direct", "time-centered", 1)

-- TODO Do automatic check.
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
