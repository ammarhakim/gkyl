-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = Moments.Eq.Euler

gasGamma = 1.4 -- gas adiabatic constant

-- create app
eulerApp = Moments.App {

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.25, 0.0}, -- lower left corner
   upper = {1.25, 2*math.pi}, -- upper right corner
   cells = {64, 64*6}, -- number of cells
   cflFrac = 0.9, -- CFL fraction

   periodicDirs = { 2 },

   mapc2p = function (t, xn)
      local r, theta = xn[1], xn[2]
      return r*math.cos(theta), r*math.sin(theta)
   end,
   
   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },

      -- initial conditions
      init = function (t, xn)
	 local r = xn[1]
	 local xs = 0.5*(0.25+1.25)

	 local rhol, ul, pl = 3, 0, 3
	 local rhor, ur, pr = 1, 0, 1

	 local rho, u, press = rhor, ur, pr
	 if r<xs then
	    rho, u, press = rhol, ul, pl
	 end
	 return rho, rho*u, 0.0, 0.0, press/(gasGamma-1) + 0.5*rho*u*u	 
      end,
      
      evolve = true, -- evolve species?

      bcx = { Moments.SpeciesBc.bcCopy, Moments.SpeciesBc.bcCopy },
   },   
}
-- run application
eulerApp:run()
