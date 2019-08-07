-- Gkyl --------------------------------------------------------------
-- BGK Relexation test -----------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

sim = Plasma.App {
   logToFile = false,

   tEnd        = 10.0,             -- End time.
   nFrame      = 1,                -- Number of frames to write.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {2*math.pi},      -- Configuration space upper right.
   cells       = {4},              -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 0.1,

   -- Decomposition for configuration space.
   decompCuts = {1},      -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},    -- Periodic directions.

   -- Neutral species (charge=1 here doesn't matter).
   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = { 6.0},
      cells = {8},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 if math.abs(v) > 1.5 then
	    return 0
	 else
	    return 1.0/3
	 end
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      -- Collisions.
      bgk = Plasma.BGKCollisions {
         collideWith     = {"neut"},
	 frequencies     = {1.0},
	 exactLagFixM012 = false,
      },
   },
}
-- Run application.
sim:run()
