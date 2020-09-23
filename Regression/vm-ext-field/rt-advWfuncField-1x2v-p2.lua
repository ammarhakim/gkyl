-- Gkyl ------------------------------------------------------------------
-- 
-- Advection in specified electromagnetic fields, non-resonant case.
-- 
--------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell()

local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd        = 50.0,          -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower left.
   upper       = {2*math.pi},   -- Configuration space upper right.
   cells       = {2},           -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4."

   -- Parallel decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- Set to true to use shared memory.

   -- Specify which directions have periodic boundary conditions for configuration space.
   periodicDirs = {1},

   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-8.0, -8.0},
      upper = {8.0, 8.0},
      cells = {16, 16},
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(1.0, vx, vy, 0.0, 0.0, 1.0)
      end,
      -- Diagnostics.
      diagnosticMoments = { "M1i" },
   },

   -- Specified, time-dependent electromagnetic fields.
   externalField = Plasma.ExternalField {
      emFunc = function (t, xn)
         -- Units for frequency and magnetic field chosen such that OmegaC = 1.0,
	 -- where OmegaC is the cyclotron frequency.
	 local omega = 0.5
	 local B0    = 1.0
	 local Ex    = 1.0*math.cos(omega*t)
	 return Ex, 0.0, 0.0, 0.0, 0.0, B0
      end,
      -- Evolve field?
      evolve = true,
   },
}
-- Run application.
sim:run()
