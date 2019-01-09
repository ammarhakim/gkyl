-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma    = require "App.PlasmaOnCartGrid"
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
me   = Constants.ELECTRON_MASS
mi   = Constants.PROTON_MASS   -- Hydrogen ions.

B0  = 0.5     -- Magnetic field amplitude [T].

n0  = 7e19    -- Number density [1/m^3].

Te0 = 40*eV   -- Electron temperature.
Ti0 = 80*eV   -- Ion temperature.

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)

-- Bulk flow speed along B-field.
uPare = math.sqrt(mi/me)*vti
uPari = 0.0

nuFrac       = 1.0
-- Electron collision frequency.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc        = nuFrac*logLambdaElc*(eV^4)*n0
              /(6*math.sqrt(2)*(math.pi^(3/2))*(eps0^2)*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision frequency.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon        = nuFrac*logLambdaIon*(eV^4)*n0
              /(12*(math.pi^(3/2))*(eps0^2)*math.sqrt(mi)*((Ti0)^(3/2)))
-- Electron-Ion collision frequency.
nuElcIon     = nuElc/1.96

-- Box size
Lx = 4 -- [m]

print('Electron-ion collision period: ', 1.0/nuElcIon)

plasmaApp = Plasma.App {
   logToFile = false,
   
   tEnd        = 1.0/nuElcIon,    -- End time.
   nFrame      = 10,              -- Number of frames to write.
   lower       = {-Lx/2},         -- Configuration space lower coordinate.
   upper       = { Lx/2},         -- Configuration space upper coordinate.
   cells       = {12},            -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 1.0,
   
   -- Decomposition for configuration space.
   decompCuts = {1},              -- Cuts in each configuration direction.
   useShared  = false,            -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},            -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   elc = Plasma.GkSpecies {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower      = {-4*vte, 0.0},
      upper      = { 6*vte, 12*me*(vte^2)/(2*B0)},
      cells      = {16, 8},
      decompCuts = {1,1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return uPare
         end,
         temperature = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return Te0
         end,
      },
      --bcx = { Plasma.GkSpecies.bcOpen,
      --        Plasma.GkSpecies.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2" },
      -- Collisions.
      coll = Plasma.GkLBOCollisions {
         collideWith = {"ion", },
         frequencies = {nuElcIon, },
      },
   },

   -- Neutral species with a bump in the tail.
   ion = Plasma.GkSpecies {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower      = {-4*vti, 0.0},
      upper      = { 4*vti, 12*mi*(vti^2)/(2*B0)},
      cells      = {16, 8},
      decompCuts = {1,1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return uPari
         end,
         temperature = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return Ti0
         end,
      },
      --bcx = { Plasma.GkSpecies.bcOpen,
      --        Plasma.GkSpecies.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2" },
      -- Collisions.
--      coll = Plasma.GkLBOCollisions {
--         collideWith = {"ion"},
--         frequencies = {nuIon},
--      },
   },

   -- Field solver.
   field = Plasma.GkField {
      evolve = true,    -- Evolve fields?
      -- initPhiFunc = function (t, xn) return 0.0 end,
      kperp2 = 0.0 
   },
   
   -- Magnetic geometry.
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },

}
-- run application
plasmaApp:run()
