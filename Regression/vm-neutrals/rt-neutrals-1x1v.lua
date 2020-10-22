-- Gkyl --------------------------------------------------------------
-- 1x1v simulation of ionization in simple test problem
-- with spatially constant fluid moments and periodic BCs.
----------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Constants = require "Lib.Constants"

-- SI units.
local eV = Constants.ELEMENTARY_CHARGE
local eps0, mu0 = Constants.EPSILON0, 1.257e-6

local qe, qi = -eV, eV
local me, mi = 9.109383e-31, 1.6726218e-27
local n0 = 1.0e17
local ne, ni, nn = n0, n0, 0.01*n0
local vde, vdi = 0.0, 0.0
local Te, Ti = 20*eV, 1*eV

local vte, vti = math.sqrt(Te/me), math.sqrt(Ti/mi)
local uB = math.sqrt(Te/mi)
local omega_pe = math.sqrt((ne * qe^2)/(eps0*me))
local lambda_D = math.sqrt((eps0 * Te)/(ne * qe^2))
-- Artificially decrease the speed of light.
mu_0 = 1.0/(eps0 * (10*vte)^2)

-- Initialization function.
local function maxwellian(n, vd, vth, v)
   return n / math.sqrt(2*math.pi*vth*vth) * 
      math.exp(-(v-vd)^2/(2*vth*vth))
end

plasmaApp = Plasma.App {
   logToFile = false,

   tEnd = 1/omega_pe,       -- End time.
   nFrame = 1,              -- Number of output frames.
   lower = {0.0},            -- Configuration space lower left.
   upper = {128.0*lambda_D}, -- Configuration space upper right.
   cells = {64},            -- Configuration space cells.
   basis = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder = 2,            -- Polynomial order.
   cflFrac = 1,              -- CFL "fraction". Usually 1.0.
   timeStepper = "rk3",      -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},         -- Cuts in each configuration direction.
   useShared = false,        -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},       -- Periodic directions.
      
   -- Electrons
   elc = Plasma.Species {
      evolve = true,
      charge = qe, mass = me,

      -- Velocity space grid.
      lower = {-6.0*vte},
      upper = {6.0*vte},
      cells = {16},
      decompCuts = {1},

      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         return maxwellian(ne, 0, vte, v)
      end,

      -- Neutral interactions
      ionization = Plasma.Ionization {
         collideWith = {"neut"},        -- species to collide with
      	 electrons = "elc",             -- define name for electron species
      	 neutrals = "neut",             -- define name for neutral species
      	 elemCharge = eV,               -- define elementary charge
      	 elcMass = me,                  -- electron mass
         plasma = "H",                  -- ion species element
      },

      -- Diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "vtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1i",
				     "intM2Flow", "intM2Thermal" },
   },

   -- Ions
   ion = Plasma.Species {
      evolve = true,
      charge = qi, mass = mi,

      -- Velocity space grid.
      lower = {-6.0*uB},
      upper = {6.0*uB},
      cells = {16},
      decompCuts = {1},

      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         return maxwellian(ni, 0, vti, v)
      end,

      -- Neutral interactions 
      ionization = Plasma.Ionization {
         collideWith = {"neut"},        -- species to collide with
      	 electrons = "elc",             -- define name for electron species
      	 neutrals = "neut",             -- define name for neutral species
      	 elemCharge = eV,               -- define elementary charge
      	 elcMass = me,                  -- electron mass
         plasma = "H",                  -- ion species element

      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"neut"},       -- species to collide with
      	 ions = "ion",                 -- define name for ion species
      	 neutrals = "neut",            -- define name for neutral species
      	 ionMass = mi,                 -- ion mass
      	 neutMass = mi,                -- neutral mass
      	 plasma = "H",                 -- ion species element
      	 charge = qi,                  -- species charge
      },

      -- Diagnostics
      diagnosticMoments = { "M0", "M1i", "M2" },
      diagnosticIntegratedMoments = {"intM0", "intM1i",
      				     "intM2Flow", "intM2Thermal" },
   },

   neut = Plasma.Species {
      evolve = true,
      charge = 0.0, mass = mi,

      -- Velocity space grid
      lower = {-6.0*uB},
      upper = {6.0*uB},
      cells = {16},
      decompCuts = {1},
      
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         return maxwellian(nn, 0, vti, v)
      end,

      -- Neutral interactions
      ionization = Plasma.Ionization {
         collideWith = {"elc"},         -- species to collide with
      	 electrons = "elc",             -- define name for electron species
      	 neutrals = "neut",             -- define name for neutral species
      	 elemCharge = eV,               -- define elementary charge
      	 elcMass = me,                  -- electron mass
         plasma = "H",                  -- ion species element  
      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"ion"},        -- species to collide with
      	 ions = "ion",                 -- define name for ion species
      	 neutrals = "neut",            -- define name for neutral species
      	 ionMass = mi,                 -- ion mass
      	 neutMass = mi,                -- neutral mass
      	 plasma = "H",                 -- ion species element
      	 charge = 0,                   -- species charge
      },

      -- Diagnostics
      diagnosticMoments = { "M0", "M1i", "M2" },
      diagnosticIntegratedMoments = {"intM0", "intM1i",
      				     "intM2Flow", "intM2Thermal" },

   },
   
   -- Field solver.
   field = Plasma.Field {
      epsilon0 = eps0, mu0 = mu0,
      init = function (t, xn)
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,
   },
}
-- Run application.
plasmaApp:run()
