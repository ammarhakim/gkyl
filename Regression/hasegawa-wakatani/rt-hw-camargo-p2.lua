-- Plasma ------------------------------------------------------------------------
--
-- Hasegawa Wakatani simulation.
--
local Plasma = require("App.PlasmaOnCartGrid").HasegawaWakatani()

local L = 2.*math.pi/0.15

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 6.0,              -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {-L/2, -L/2},     -- Configuration space lower left.
   upper       = { L/2,  L/2},     -- Configuration space upper right.
   cells       = {172,172},        -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 1.,
   
   -- Decomposition for configuration space.
   decompCuts = {1,1},    -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2},    -- Periodic directions.

   -- Fluid species.
   fluid = Plasma.Species {
      adiabaticity = 1.0,  gradient = 1.0,
      -- Initial conditions.
      randomseed = 1234,  -- So regression test doesn't fail.
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local r0 = 0.1
         local vorticity, density = 0., 2*r0*(math.random()-0.5)
         return vorticity, density
      end,
      diff = Plasma.Diffusion {
         order = 6,
         coefficient = 1.e-8,
      },
      evolve = true,    -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true,    -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
