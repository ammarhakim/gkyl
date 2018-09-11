-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

app = Plasma.App {
   logToFile = false,

   tEnd = 1.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {0}, -- configuration space lower left
   upper = {1}, -- configuration space upper right
   cells = {2}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   noFix = Plasma.VlasovMaxwell.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0,},
      upper = {6.0,},
      cells = {4,},
      decompCuts = {1,},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = 1.0,
         drift = 1.0,
         temperature = 1.0,
         exactScaleM0 = false,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   scaleFix = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0,},
      upper = {6.0,},
      cells = {4,},
      decompCuts = {1,},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = 1.0,
         drift = 1.0,
         temperature = 1.0,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   lagFix = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0,},
      upper = {6.0,},
      cells = {4,},
      decompCuts = {1,},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = 1.0,
         drift = 1.0,
         temperature = 1.0,
         exactScaleM0 = false,
         exactLagFixM012 = true,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   double = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0,},
      upper = {6.0,},
      cells = {4,},
      decompCuts = {1,},
      -- initial conditions
      init1 = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = 1.0,
         drift = 1.0,
         temperature = 1.0,
         exactScaleM0 = false,
         exactLagFixM012 = true,
      },
      init2 = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = 1.0,
         drift = -1.0,
         temperature = 1.0,
         exactScaleM0 = false,
         exactLagFixM012 = true,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },
}
