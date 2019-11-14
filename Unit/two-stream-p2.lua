-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- This is the same two-stream p2 simulation as in Regression.
knumber      = 0.5    -- Wave-number.
elVTerm      = 0.2    -- Electron thermal velocity.
vDrift       = 1.0    -- Drift velocity.
perturbation = 1.0e-6 -- Distribution function perturbation.
elVTerm      = 0.2    -- Electron thermal velocity.

vlasovApp = Plasma.App {
   logToFile = false,

   tEnd        = endTime,            -- End time.
   nFrame      = frames,             -- Number of output frames.
   lower       = {-math.pi/knumber}, -- Configuration space lower left.
   upper       = { math.pi/knumber}, -- Configuration space upper right.
   cells       = {64},               -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = decomp,   -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space
   periodicDirs = {1}, -- Periodic directions.

   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = {6.0},
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
         local x, v  = xn[1], xn[2]
         local alpha = perturbation
         local k     = knumber
         local vt    = elVTerm

         local fv = 1/math.sqrt(8*math.pi*vt^2)*(math.exp(-(v-vDrift )^2/(2*vt^2))+math.exp(-(v+vDrift)^2/(2*vt^2)))
         return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- Evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k     = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
