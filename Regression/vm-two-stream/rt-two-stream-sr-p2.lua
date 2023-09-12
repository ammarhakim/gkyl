-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

knumber      = 0.5    -- Wave-number.
temp         = 0.04   -- Electron temperature
vDrift       = 0.9    -- Drift velocity.
perturbation = 1.0e-3 -- Distribution function perturbation.

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd        = 100.0,              -- End time.
   nFrame      = 1,                  -- Number of output frames.
   lower       = {-math.pi/knumber}, -- Configuration space lower left.
   upper       = {math.pi/knumber},  -- Configuration space upper right.
   cells       = {64},              -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space
   periodicDirs = {1}, -- Periodic directions.

   -- Electrons.
   elc = Plasma.GenSpecies.VlasovSR {
      charge = -1.0,  mass = 1.0,
      -- Velocity space grid.
      lower = {-8.0},
      upper = { 8.0},
      cells = {64},
      -- Initial conditions.
      init = function (t, xn)
         local x, p = xn[1], xn[2]

         local K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12
         local gamma = 1.0/math.sqrt(1 - vDrift*vDrift)
         local alpha = perturbation
         local k = knumber
         local n = 1+alpha*math.cos(k*x)
         local mc2_T = 1.0/temp

         local fv = 0.5*n/(4*math.pi*temp*K_2)*(math.exp(-mc2_T*gamma*(math.sqrt(1 + p*p) - vDrift*p)) + math.exp(-mc2_T*gamma*(math.sqrt(1 + p*p) + vDrift*p)))
         return fv
      end,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i" }
   },

   -- Field solver
   field = Plasma.GenField.Maxwell {
      epsilon0 = 1.0,  mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
