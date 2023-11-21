-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic()
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
c_s = math.sqrt(Te0/mi)
omega_ci = math.sqrt(qi*B0/mi)
rho_s = c_s / omega_ci
--print(string.format("c_s=%10.8e, omega_ci=%10.8e, rho_s=%10.8e\n", c_s, omega_ci, rho_s))

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
--print(string.format("nuElc=%10.8e, nuIon=%10.8e, nuElcIon=%10.8e, nuIonElc=%10.8e", nuElc, nuIon, nuElcIon, nuIonElc))

-- Box size.
Lx = 4 -- [m]
Ly = 4
Lz = 4 -- [m]

plasmaApp = Plasma.App {
   logToFile = false,
   
   tEnd        = 0.10/nuElcIon,   -- End time.
   nFrame      = 1,               -- Number of frames to write.
   lower       = {-Lx/2, -Ly/2, -Lz/2},   -- Configuration space lower coordinate.
   upper       = { Lx/2,  Ly/2,  Lz/2},   -- Configuration space upper coordinate.
   cells       = {4, 4, 4},       -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 1.0,
   
   -- Decomposition for configuration space.
   decompCuts = {4,4,4},              -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1,2,3},            -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower      = {-5*vte, 0.0},
      upper      = { 5*vte, 12*me*(vte^2)/(2*B0)},
      cells      = {8, 4},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return n0 * (1.0+0.2*math.cos(2.*math.pi*x))
         end,
         driftSpeed = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return uPare * (1.0+0.2*math.cos(2.*math.pi*y))
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return Te0 * (1.0+0.2*math.cos(2.*math.pi*z))
         end,
      },
      polarizationDensityFactor = n0,
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnostics = { "M0", "M1", "M2", "Upar", "VtSq", "intM0", "intM1", "intM2" },
      -- Collisions.
      coll = Plasma.BGKCollisions {
         collideWith = { "elc", "ion" },
         frequencies = { nuElc, nuElcIon },
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
      cells      = {8, 4},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return n0 * (1.0+0.2*math.cos(2.*math.pi*x))
         end,
         driftSpeed = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return uPari * (1.0+0.2*math.cos(2.*math.pi*y))
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return Ti0 * (1.0+0.2*math.cos(2.*math.pi*z))
         end,
      },
      polarizationDensityFactor = n0,
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnostics = { "M0", "M1", "M2", "Upar", "VtSq", "intM0", "intM1", "intM2" },
      -- Collisions.
      coll = Plasma.BGKCollisions {
         collideWith = { "ion", "elc" },
         frequencies = { nuIon, nuIonElc },
         -- Optional arguments:
--         betaGreene  = 1.0,    -- Free parameter, must be >-1.
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = false,    -- Evolve fields?
      externalPhi = function (t, xn) return 0.0 end,
      kperpSq = 0.0 
   },
   
   -- Magnetic geometry.
   externalField = Plasma.Geometry {
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
