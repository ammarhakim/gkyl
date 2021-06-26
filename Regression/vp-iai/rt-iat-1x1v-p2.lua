-- Gkyl -----------------------------------------------------------------------
--
-- Ion acoustic instability simulation.
--
-------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Electron parameters.
vDriftElc = 0.01   -- Modified from 0.159.
vtElc     = 0.02
-- Ion parameters.
vDriftIon = 0.0
vtIon     = 0.000566      -- Modified from 0.001 (use Te = 50 Ti).
-- Mass ratio.
massRatio = 25.0  -- Modified from 25.

knumber      = 10.0    -- Wave-number.
perturbation = 1.0e-4  -- Distribution function perturbation.

local function maxwellian1v(v, vDrift, vt)
   return 1/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vDrift)^2/(2*vt^2))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 30.0,          -- End time.
   nFrame      = 1,             -- Number of output frames.

   lower       = {0.0},         -- Configuration space lower left.
   upper       = {1.0},         -- Configuration space upper right.
   cells       = {64},          -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2" or "rk3".
   cflFrac     = 0.9,

   -- Eecomposition for configuration space
   decompCuts = {1},  -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.
   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.001,   

   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- Velocity space grid.
      lower      = {-6.0*vDriftElc},
      upper      = { 6.0*vDriftElc},
      cells      = {128},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         local fv = maxwellian1v(v, vDriftElc, vtElc)
	 return fv*(1)
      end,
      evolve = true, -- Evolve species?

      diagnostics = { "M0", "M1i", "intM0", "intM1i" },
   },

   -- Ions.
   ion = Plasma.Species {
      charge = 1.0, mass = massRatio,
      -- Velocity space grid.
      lower      = {-6.0*vtIon},   
      upper      = { 6.0*vtIon},
      cells      = {128},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         local fv = maxwellian1v(v, vDriftIon, vtIon)
         return fv*(1+perturbation*math.cos(2*math.pi*knumber*x))
      end,
      evolve = true, -- Evolve species?

      diagnostics = { "M0", "M1i", "intM0", "intM1i" },
   },   

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0,
      evolve   = true, -- Evolve field?
      hasMagneticField = false,
   },

}
-- Run application.
plasmaApp:run()
