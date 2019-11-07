-- Plasma ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").IncompEuler

local L = 10.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 100,            -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0, 0},           -- Configuration space lower left.
   upper       = {L, L},           -- Configuration space upper right.
   cells       = {64,64},          -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 1.,
   
   -- Decomposition for configuration space.
   decompCuts = {1,1},    -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2}, -- Periodic directions.

   -- Fluid species.
   fluid = Plasma.Species {
      charge = -1.0,
      -- Initial conditions.
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local x1, y1 = 3.5, 5.0
         local x2, y2 = 6.5, 5.0
         local r1 = (x-x1)^2 + (y-y1)^2
         local r2 = (x-x2)^2 + (y-y2)^2
         return math.exp(-r1/0.8^2) + math.exp(-r2/0.8^2) -- 4.0212385953656/L^2
      end,
      evolve          = true, -- Evolve species?
      applyPositivity = true,
      positivityDiffuse = true,
      positivityRescale = true,
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
