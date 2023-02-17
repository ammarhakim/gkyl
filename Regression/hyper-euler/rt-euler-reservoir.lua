-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()

gasGamma = 3.0 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {
   tEnd = 50.0, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {2}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   
   -- decomposition stuff
   decompCuts = {1}, -- cuts in each direction

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 return 1.0, 0.0, 0.0, 0.0, 0.1
      end,
      
      evolve = true, -- evolve species?

      bcx = { Moments.Species.bcConst(1.0, 0.0, 0.0, 0.0, 1.0/(gasGamma-1)),
	      Moments.Species.bcConst(0.125, 0.0, 0.0, 0.0, 0.1/(gasGamma-1)) },
	      
   },   
}
-- run application
eulerApp:run()
