-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local Constants = require "Lib.Constants"
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = True
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

----------------
-- PARAMETERS --
----------------
-- physical constants and normalizations
             lightSpeed = Constants.SPEED_OF_LIGHT
                    mu0 = Constants.MU0
               gasGamma = 5. / 3.
                     me = Constants.ELECTRON_MASS
                     qe = -Constants.ELEMENTARY_CHARGE
                     mi = 131.293 * Constants.PROTON_MASS
                     qi = Constants.ELEMENTARY_CHARGE
                     kB = Constants.BOLTZMANN_CONSTANT

               epsilon0 = 1. / mu0 / lightSpeed^2
-- problem setup
                     E0 = 20e3  -- along y
                     B0 = 0.02  -- along z
                     n0 = 1e17
                    Te0 = 10 * Constants.EV2KELVIN
                    Ti0 = 0.2 * Constants.EV2KELVIN
                      L = 44.56e-3
                   pert = 1e-7  -- perturbation level

                    Te0 = 10 * Constants.EV2KELVIN * 0.25
                      L = 2.82e-3

                  vExB0 = E0 / B0  -- along x
                    pe0 = n0 * kB * Te0
                    pi0 = n0 * kB * Ti0
-- diagnostic parameters for verification
                   vte0 = math.sqrt(kB * Te0 / me)
                   vti0 = math.sqrt(kB * Ti0 / mi)
              lambdaDe0 = math.sqrt(epsilon0 * kB * Te0 / n0 / qe^2)
              lambdaDi0 = math.sqrt(epsilon0 * kB * Ti0 / n0 / qi^2)
               lambdaD0 = math.sqrt(epsilon0 * kB / n0 / qe^2 / (1. / Te0 + 1. / Ti0))
                   wpe0 = math.sqrt(n0 * qe^2 / me / epsilon0)
                   wce0 = math.abs(qe * B0 / me)
                   wpi0 = math.sqrt(n0 * qi^2 / mi / epsilon0)
                   wci0 = math.abs(qi * B0 / mi)
                    de0 = lightSpeed / wpe0
                    di0 = lightSpeed / wpi0
                  pmag0 = B0^2 / 2. / mu0
                     p0 = (pe0 + pi0)
                  beta0 = p0/ pmag0
-- domain and grid
                     Lx = L
                     Nx = 256
                  lower = {0.}
                  upper = {Lx}
                  cells = {Nx}
           periodicDirs = {1}
-- computational parameters
              numFluids = 3
                 charge = {qe, qi, qe}
                   mass = {me, mi, me}
                    cfl = 0.9
                limiter = "monotonized-centered"
-- i/o control
                   tEnd = 1000e-9
                 tFrame = 1e-9
                 nFrame = tEnd / tFrame

log("%30s = %g", "mi/me", mi / me)
log("%30s = %g", "wpe0/wce0", wpe0 / wce0)
log("%30s = %g", "beta0", beta0)
log("%30s = %g", "pi0/pe0", pi0 / pe0)
log("%30s = %g", "vExB0/c", vExB0 / lightSpeed)
log("%30s = %g", "vte0/c", vte0 / lightSpeed)
log("%30s = %g", "vti0/c", vti0 / lightSpeed)
log("%30s = %g", "wpe0", wpe0)
log("%30s = %g", "wce0", wce0)
log("%30s = %g", "wpi0", wpi0)
log("%30s = %g", "wci0", wci0)
log("%30s = %g", "de0", de0)
log("%30s = %g", "di0", di0)
log("%30s = %g", "de0/L", de0 / L)
log("%30s = %g", "di0/L", di0 / L)
log("%30s = %g", "lambdaDe0", lambdaDe0)
log("%30s = %g", "lambdaDi0", lambdaDi0)
log("%30s = %g", "lambdaD0", lambdaD0)
log("%30s = %g", "lambdaDe0/L", lambdaDe0 / L)
log("%30s = %g", "lambdaDi0/L", lambdaDi0 / L)
log("%30s = %g", "lambdaD0/L", lambdaD0 / L)
log("%30s = %g", "pert", pert)

-----------------------
-- INITIAL CONDITION --
-----------------------
math.randomseed(os.time())
m_max = math.floor(Nx / 4)
m_max = 36
waveModes = {}
for m = 1, m_max do
  waveModes[m] = m
end
waveNumbers = {}
for i, m in ipairs(waveModes) do
   waveNumbers[i] = 2. * math.pi / L * m
end

wavePhases = {}
math.randomseed(os.time())
for i, m in ipairs(waveModes) do -- random phase for each wave
   wavePhases[i] = 2. * math.pi * m * math.random()
end

init = function(x, y, z)
   local n = n0
   local vx_e, vy_e, vz_e = vExB0, 0., 0.
   local vx_i, vy_i, vz_i = 0., 0., 0.
   local p_e = pe0
   local p_i = pi0
   local Ex, Ey, Ez = 0., E0, 0.
   local Bx, By, Bz = 0., 0., B0
   
   local n_e = n
   local n_i = n
   for m,k_m in ipairs(waveNumbers) do
      -- local dn_e = n * pert
      -- local dEx =  - qe * dn_e / k_m / epsilon0
      local dn_e1 = n * pert
      local dEx =  - qe * dn_e1 / epsilon0
      local dn_e = dn_e1 * k_m
      -- local dn_e = - dEx * k_m * epsilon0 / qe
      n_e = n_e + dn_e * math.sin(k_m * x + wavePhases[m])
      Ex = Ex + dEx * math.cos(k_m * x + wavePhases[m])
   end
   
   local rho_e = n_e * me
   local rho_i = n_i * mi
   local rhovx_e = rho_e * vx_e
   local rhovy_e = rho_e * vy_e
   local rhovz_e = rho_e * vz_e
   local rhovx_i = rho_i * vx_i
   local rhovy_i = rho_i * vy_i
   local rhovz_i = rho_i * vz_i
   local e_e = p_e / (gasGamma - 1.) + 0.5 * (rhovx_e^2 + rhovy_e^2 + rhovz_e^2) / rho_e
   local e_i = p_i / (gasGamma - 1.) + 0.5 * (rhovx_i^2 + rhovy_i^2 + rhovz_i^2) / rho_i
   
   -- ghost species
   local rho_g = n * me
   local rhovx_g = - rho_g * vx_e
   local rhovy_g = - rho_g * vy_e
   local rhovz_g = - rho_g * vz_e
   local e_g = e_e

   return rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
          rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
          rho_g, rhovx_g, rhovy_g, rhovz_g, e_g,
          Ex, Ey, Ez, Bx, By, Bz, 0., 0.
end

----------------
-- GKEYLL App --
----------------
momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = lower,
   upper = upper,
   cells = cells,
   timeStepper = "fvDimSplit",

   decompCuts = {4}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   periodicDirs = periodicDirs, -- periodic directions

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]

         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
            rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
            rho_g, rhovx_g, rhovy_g, rhovz_g, e_g,
            Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x)
      
          return rho_e, rhovx_e, rhovy_e, rhovz_e, e_e
      end,
      evolve = true, -- evolve species?
      limiter = limiter,
      cfl = cfl,
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]
         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
            rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
            rho_g, rhovx_g, rhovy_g, rhovz_g, e_g,
            Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x)
      
          return rho_i, rhovx_i, rhovy_i, rhovz_i, e_i
      end,
      evolve = true, -- evolve species?
      limiter = limiter,
      cfl = cfl,
   },

   -- ghost electrons
   gho = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]

         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
            rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
            rho_g, rhovx_g, rhovy_g, rhovz_g, e_g,
            Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x)
      
          return rho_g, rhovx_g, rhovy_g, rhovz_g, e_g
      end,
      evolve = false, -- evolve species?
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x = xn[1]

         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
            rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
            rho_g, rhovx_g, rhovy_g, rhovz_g, e_g,
            Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x)
      
          return Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
      end,
      evolve = true, -- evolve field?
      limiter = limiter,
      cfl = cfl,
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion", "gho"},
      timeStepper = "analytic2",
   },   

}
-- run application
momentApp:run()
