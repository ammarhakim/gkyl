-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local SrEuler = Moments.Eq.SrEuler

gasGamma = 5.0/3.0 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {

   tEnd = 50.0, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {100.0}, -- upper right corner
   cells = {1024}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   
   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = SrEuler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local xs = 45.0

	 local rhol, ul, pl = 10.0, 0.0, 40./3.;
	 local rhor, ur, pr = 1.0, 0.0, 2./(3.e7);

	 local rho, u, p = rhor, ur, pr
	 if xn[1]<xs then
	    rho, u, p = rhol, ul, pl
	 end

	 local gamma = 1 / math.sqrt(1 - u*u);
	 local rhoh = gasGamma * p / (gasGamma - 1)  + rho;

	 return gamma*rho, gamma*gamma*rhoh - p, gamma*gamma*rhoh*u, 0.0, 0.0
      end,
      evolve = true, -- evolve species?

      bcx = { Moments.SpeciesBc.bcCopy, Moments.SpeciesBc.bcCopy }
   },   
}
-- run application
eulerApp:run()
