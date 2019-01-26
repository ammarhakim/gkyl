-- Gkyl ------------------------------------------------------------------------
--
-- Test for five-moment source updater
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local epsilon0 = 1.0
local charge = {}
local mass ={}
local rho = {}
local rhovx = {}
local rhovy = {}
local rhovz = {}
local pressure = {}

charge[1], mass[1], rho1, rhovx1, rhovy1, rhovz1, p1 = -1, 1, 4, 1, 2, 3, 1
charge[2], mass[2], rho2, rhovx2, rhovy2, rhovz2, p2 = 1, 4, 16, 6, 8, 10, 1.5

rho = {rho1, rho2}

local Ex, Ey, Ez = 1.2, 1.3, 1.4

-- local Bx, By, Bz = 0, 0, 0  -- zero B field, parallel sol. only
-- timesteps = {0.01}

-- local Bx, By, Bz = 0, 0, 3.
-- timesteps = {0.01}

-- local Bx, By, Bz = 0, 3, 0
-- timesteps = {0.01}

-- local Bx, By, Bz = 3, 0, 0
-- timesteps = {0.01}

local Bx, By, Bz = 1, 2, 3

local Bmag = math.sqrt(Bx^2 + By^2 + Bz^2)

for i = 1,#charge do
   print(string.format("[%d] charge=%g, mass=%g, rho=%g",
      i, charge[i], mass[i], rho[i]))
end

wp_tot = 0.
wc_max = 0.
for i,_ in ipairs(charge)  do
   local wp2 = rho[i] * (charge[i] / mass[i])^2 / epsilon0
   wp_tot = wp_tot + wp2
   wc_max = math.max(wc_max, math.abs(charge[i] * Bmag / mass[i]))
end
wp_tot = math.sqrt(wp_tot)
print('wp_tot', wp_tot)
print('plasma period', 2 * math.pi / wp_tot)
print('wc_max', wc_max)
print('cyclotron period', 2 * math.pi / wc_max)

timesteps = {0.001*math.pi/wp_tot, 200.001*math.pi/wp_tot}

local ffi  = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

local grid = Grid.RectCart {
   lower = {0.0},
   upper = {1.0},
   cells = {1},
}
local fluid1 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   ghost = {1, 1},
}
local fluid2 = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   ghost = {1, 1},
}
local em = DataStruct.Field {
   onGrid = grid,
   numComponents = 8,
   ghost = {1, 1}
}

local localRange = fluid1:localRange()
local fluid1Idxr = fluid1:indexer()   
local fluid2Idxr = fluid2:indexer()   
local emIdxr = em:indexer() 

local gasGamma = 5. / 3.

function init()
   for i = localRange:lower(1), localRange:upper(1) do
      local fluid1Itr = fluid1:get(fluid1Idxr(i))
      local fluid2Itr = fluid2:get(fluid2Idxr(i))
      local emItr = em:get(emIdxr(i))

      local e1 = 0.5 * (rhovx1^2 + rhovy1^2 + rhovz1^2) / rho1 + p1 / (gasGamma - 1)
      fluid1Itr[1] = rho1
      fluid1Itr[2] = rhovx1
      fluid1Itr[3] = rhovy1
      fluid1Itr[4] = rhovz1
      fluid1Itr[5] = e1

      local e2 = 0.5 * (rhovx2^2 + rhovy2^2 + rhovz2^2) / rho2 + p2 / (gasGamma - 1)
      fluid2Itr[1] = rho2
      fluid2Itr[2] = rhovx2
      fluid2Itr[3] = rhovy2
      fluid2Itr[4] = rhovz2
      fluid2Itr[5] = e2

      emItr[1] = Ex
      emItr[2] = Ey
      emItr[3] = Ez
      emItr[4] = Bx
      emItr[5] = By
      emItr[6] = Bz

      emItr[7] = 3.33
      emItr[8] = 4.56
   end
end

function display()
   for i = localRange:lower(1), localRange:upper(1) do
      local fluid1Itr = fluid1:get(fluid1Idxr(i))
      local fluid2Itr = fluid2:get(fluid1Idxr(i))
      local emItr = em:get(emIdxr(i))
      
      print(string.format("fluid1 = {%.8e, %.8e, %.8e, %.8e, %.8e}", fluid1Itr[1], fluid1Itr[2], fluid1Itr[3], fluid1Itr[4], fluid1Itr[5]))
      print(string.format("fluid2 = {%.8e, %.8e, %.8e, %.8e, %.8e}", fluid2Itr[1], fluid2Itr[2], fluid2Itr[3], fluid2Itr[4], fluid2Itr[5]))
      print(string.format("ele = {%.8e, %.8e, %.8e}", emItr[1], emItr[2], emItr[3]))
      print(string.format("mag = {%g, %g, %g}", emItr[4], emItr[5], emItr[6]))
      print(string.format("phiE, phiM = {%g, %g}", emItr[7], emItr[8]))
   end
   print("\n")
end

function test_1(scheme, dt)
   init()
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
   srcUpdater:advance(0.0, {}, {fluid1, fluid2, em})
   display()
end

init()
print("----- BEFORE -----")
display()

for _,scheme in ipairs({"analytic2", "exact"}) do
   for _,dt in ipairs(timesteps) do
      print(string.format("+++++ AFTER; scheme = %12s, dt = %g = %g pi/wp_tot +++++",
                         scheme, dt, dt*wp_tot/math.pi))
      test_1(scheme, dt)
   end
end

-- TODO Do automatic check.
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
