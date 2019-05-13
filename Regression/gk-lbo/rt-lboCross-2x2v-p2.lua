-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic
local Constants = require "Lib.Constants"

-- This test initializes Maxwellian electrons and ions with different
-- bulk velocity and temperature and collides them.
-- We will use cross collisions only as a default, but one can add
-- self collisions as well.

-- Universal parameters.
eps0 = Constants.EPSILON0
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV
qi   =  eV
me   = Constants.ELECTRON_MASS --*1836.0
mi   = Constants.PROTON_MASS   -- Hydrogen ions.

B0  = 0.5     -- Magnetic field amplitude [T].

n0  = 7e19    -- Number density [1/m^3].

Te0 = 300*eV   -- Electron temperature.
Ti0 = 200*eV   -- Ion temperature.

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)

-- Bulk flow speed along B-field.
uPare = 0.5*math.sqrt(me/mi)*vte
uPari = 50.0*(me/mi)*vti

nuFrac       = 1.0
-- Electron collision frequency.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc        = nuFrac*logLambdaElc*(eV^4)*n0
              /(12*math.sqrt(2)*(math.pi^(3/2))*(eps0^2)*math.sqrt(me)*(Te0^(3/2)))
-- Ion collision frequency.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon        = nuFrac*logLambdaIon*(eV^4)*n0
              /(12*math.sqrt(2)*(math.pi^(3/2))*(eps0^2)*math.sqrt(mi)*(Ti0^(3/2)))
-- Electron-Ion collision frequency.
nuElcIon     = nuElc*2.0
nuIonElc     = me*nuElcIon/mi -- Ion-electron collision frequency.

-- Box size.
Lx = 4 -- [m]
Ly = 4 -- [m]

plasmaApp = Plasma.App {
   logToFile = false,
   
   tEnd        = 0.01/nuElcIon,   -- End time.
   nFrame      = 1,               -- Number of frames to write.
   lower       = {-Lx/2, -Ly/2},  -- Configuration space lower coordinate.
   upper       = { Lx/2,  Ly/2},  -- Configuration space upper coordinate.
   cells       = {12, 12},        -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 1.0,
   
   -- Decomposition for configuration space.
   decompCuts = {1, 1},           -- Cuts in each configuration direction.
   useShared  = false,            -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2},         -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower      = {-5*vte, 0.0},
      upper      = { 5*vte, 12*me*(vte^2)/(2*B0)},
      cells      = {16, 8},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return uPare
         end,
         temperature = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return Te0
         end,
      },
      --bcx = { Plasma.GkSpecies.bcOpen,
      --        Plasma.GkSpecies.bcOpen },
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1", "intM2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = { "elc", "ion" },
         frequencies = { nuElc, nuElcIon },
--         collideWith = { "elc" },
--         frequencies = { nuElc },
--         collideWith = { "ion" },
--         frequencies = { nuElcIon },
         -- Optional arguments:
--         betaGreene  = 1.0,    -- Free parameter, must be >-1.
      },
   },

   -- Neutral species with a bump in the tail.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower      = {-5*vti, 0.0},
      upper      = { 5*vti, 12*mi*(vti^2)/(2*B0)},
      cells      = {16, 8},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return uPari
         end,
         temperature = function (t, xn)
            local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
            return Ti0
         end,
      },
      --bcx = { Plasma.GkSpecies.bcOpen,
      --        Plasma.GkSpecies.bcOpen },
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1", "intM2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = { "ion", "elc" },
         frequencies = { nuIon, nuIonElc },
--         collideWith = { "ion" },
--         frequencies = { nuIon },
--         collideWith = { "elc" },
--         frequencies = { nuIonElc },
         -- Optional arguments:
--         betaGreene  = 1.0,    -- Free parameter, must be >-1.
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = false,    -- Evolve fields?
--      initPhiFunc = function (t, xn) return 0.0 end,
      kperp2 = 0.0,
   },
   
   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,
      -- Geometry is not time-dependent.
      evolve = false,
   },

}
-- Run application.
plasmaApp:run()
