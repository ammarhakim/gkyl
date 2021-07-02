-- Gkyl ------------------------------------------------------------------------
--
-- Apply diffusion to a 1D function, using the incompressible Euler solver.
--
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").IncompEuler()

local diffCoefficient = 0.01

-- Options to benchmark against Van Leer (2005) using projection.
local function ICfunc(t, xn)
   local x = xn[1]
   local k = 2.0*math.pi
   return math.cos(k*x)
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 4.0,            -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {1.0},            -- Configuration space upper right.
   cells       = {16},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3opSplit",
--   timeStepper = "rk3",
--   maximumDt = 0.1,
--   suggestedDt = 0.1,
   
   periodicDirs = {1},

   -- Decomposition for configuration space.
   decompCuts = {1},      -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Fluid species.
   fluid = Plasma.Species {
      charge = 1.0,
      -- Projected initial conditions (with quadrature).
      init   = ICfunc,
      evolve              = true, -- Evolve species?
      evolveCollisionless = false,
      diff = Plasma.Diffusion {
         coefficient = diffCoefficient,
         treatment = "sts",
      },
   },

}
-- Run application.
plasmaApp:run()
