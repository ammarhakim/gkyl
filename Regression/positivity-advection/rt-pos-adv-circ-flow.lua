-- Plasma ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").IncompEuler

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 4*math.pi,        -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0, 0},           -- Configuration space lower left.
   upper       = {1.0, 1.0},       -- Configuration space upper right.
   cells       = {32, 32},         -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   
   -- Decomposition for configuration space.
   decompCuts = {1, 1},    -- Cuts in each configuration direction.
   useShared  = false,     -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2},    -- Periodic directions.

   -- Fluid species.
   fluid = Plasma.Species {
      charge = 1.0,
      -- Initial conditions.
      init = function (t, xn)
	 local x, y       = xn[1], xn[2]
	 local x0, y0, r0 = 0.25, 0.5, 0.15
	 local r          = math.min(math.sqrt((x-x0)^2+(y-y0)^2), r0)/r0
	 return 0.25*(1+math.cos(math.pi*r))
      end,
      evolve          = true, -- Evolve species?
      applyPositivity = true,
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = false,    -- Evolve field?
      -- u = {dphi/dy, -dphi/dx}
      initPhiFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         return -0.5*(y^2-y+x^2-x)
      end, 
   },
}
-- Run application.
plasmaApp:run()
