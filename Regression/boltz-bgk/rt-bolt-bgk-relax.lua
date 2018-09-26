-- Gkyl --------------------------------------------------------------
-- BGK Relexation test -----------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

sim = Plasma.App {
   logToFile = false,

   tEnd = 10.0, -- end time
   nFrame = 1, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {2*math.pi}, -- configuration space upper right
   cells = {4}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"
   cflFrac = 0.1,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {8},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 if math.abs(v) > 1.5 then
	    return 0
	 else
	    return 1.0/3
	 end
      end,
       -- evolve species?
      evolve = true,
      -- diagnostic moments
      diagnosticMoments = { "M0", "M1i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i",
				      "intM2Flow", "intM2Thermal" },
         -- collisions
      bgk = Plasma.BgkCollisions {
	 collFreq = 1.0,
      },
   },
}
-- run application
sim:run()
