-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 3.0 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 50.0, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {2}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 return 1.0, 0.0, 0.0, 0.0, 0.1
      end,
      
      evolve = true, -- evolve species?

      bcx = { Euler.bcConst(1.0, 0.0, 0.0, 0.0, 1.0/(gasGamma-1)),
	      Euler.bcConst(0.125, 0.0, 0.0, 0.0, 0.1/(gasGamma-1)) },
	      
   },   
}
-- run application
eulerApp:run()
