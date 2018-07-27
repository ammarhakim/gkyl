-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Euler = require "Eq.Euler"

gasGamma = 1.4 -- gas adiabatic constant
   
-- create app
eulerApp = Plasma.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {1024}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- electrons
   fluid = Plasma.Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local xs = 0.5
	 local rhol, ul, pl = 3, 0, 3
	 local rhor, ur, pr = 1, 0, 1

	 local rho, u, press = rhor, ur, pr
	 if xn[1]<xs then
	    rho, u, press = rhol, ul, pl
	 end
	 return rho, rho*u, 0.0, 0.0, press/(gasGamma-1) + 0.5*rho*u*u
      end,
      evolve = true, -- evolve species?
   },   
}
-- run application
eulerApp:run()
