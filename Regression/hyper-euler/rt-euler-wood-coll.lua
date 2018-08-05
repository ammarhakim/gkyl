-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 1.4 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.038, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {400}, -- number of cells
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
	 local x = xn[1]
	 local x1, x2 = 0.1, 0.9 -- location of discontinuities
	 local pl, pm, pr = 1000.0, 0.01, 100.0 -- pressure in each region
	 local press = 0.0
	 if x<x1 then
	    press = pl
	 elseif x<x2 then
	    press = pm
	 else
	    press = pr
	 end
	 local er = press/(gasGamma-1)
	 return 1.0, 0.0, 0.0, 0.0, er
      end,
      evolve = true, -- evolve species?

      -- boundary conditions in X
      bcx = { Euler.bcWall, Euler.bcWall },
   },   
}
-- run application
eulerApp:run()
