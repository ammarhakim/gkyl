-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------

local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic
local Constants = require "Lib.Constants"

-- Universal parameters.
eps0 = Constants.EPSILON0
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV
qi   =  eV
me   = Constants.ELECTRON_MASS
mi   = Constants.PROTON_MASS   -- Hydrogen ions.

B0  = 0.5     -- Magnetic field amplitude [T].

n0  = 1e17    -- Number density [1/m^3].

Te0 = 20*eV   -- Electron temperature.
Ti0 = 20*eV   -- Ion temperature.
Tn0 = 1*eV

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)
vtn = math.sqrt(Tn0/mi)
cs = math.sqrt(Te0/mi)

Lx = 4 --[m]

sim = Plasma.App {
   logToFile = false,

   tEnd        = 1*Lx/vte,  --1000/omega_pe,    -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {-Lx/2},            -- Configuration space lower left.
   upper       = {Lx/2}, -- Configuration space upper right.
   cells       = {32},            -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   cflFrac     = 1,                -- CFL "fraction". Usually 1.0.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.
      
    -- Electrons.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {-6.0*vte},
      upper = {6.0*vte},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 0 
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Te0
         end,
      },
      evolve = true,
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
				     "intM2" },
   },

      -- Ions
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {-6.0*vti},
      upper = {6.0*vti},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return cs
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Ti0
         end,
      },
      evolve = true,
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
				     "intM2" },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"neut"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      }
   },

   neut = Plasma.Species {
      charge = 0.0, mass = mi,
      -- Velocity space grid
      lower = {-6.0*vti},
      upper = {6.0*vti},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = {"maxwellian",
         density = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 100*n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 0
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Ti0
         end,
      },
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
				     "intM2" },
      evolve = true, -- Evolve species?
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"ion"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      }
   },
   
    -- Field solver.
   field = Plasma.Field {
      evolve = true,
      kperp2 = 0.0
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- Run application.
sim:run()
