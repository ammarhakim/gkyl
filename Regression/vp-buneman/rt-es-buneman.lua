-- Gkyl ------------------------------------------------------------------------
--
-- A 1x1v Buneman instability simulation.
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Electron parameters.
vDriftElc = 0.159
vtElc     = 0.02
-- Ion parameters.
vDriftIon = 0.0
vtIon     = 0.001
-- Mass ratio.
massRatio = 25.0

knumber      = 1.0    -- Wave-number.
perturbation = 1.0e-6 -- Distribution function perturbation.

local function maxwellian1v(v, vDrift, vt)
   return 1/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vDrift)^2/(2*vt^2))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 60.0,          -- End time.
   nFrame      = 1,             -- Number of output frames.
   lower       = {0.0},         -- Configuration space lower left.
   upper       = {1.0},         -- Configuration space upper right.
   cells       = {16},          -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions

   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0,  mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0*vDriftElc},
      upper = { 6.0*vDriftElc},
      cells = {64},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local fv = maxwellian1v(v, vDriftElc, vtElc)
	 return fv*(1+perturbation*math.cos(2*math.pi*knumber*x))
      end,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i", "M2" }
   },

   -- Electrons
   ion = Plasma.Species {
      charge = 1.0,  mass = massRatio,
      -- Velocity space grid.
      lower = {-32.0*vtIon},
      upper = { 32.0*vtIon},
      cells = {64},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian1v(v, vDriftIon, vtIon)
      end,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i", "M2" }
   },   

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0,  mu0 = 1.0,
      hasMagneticField = false,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
