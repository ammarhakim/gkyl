-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = Moments.Eq.Euler

gasGamma = 1.4 -- gas adiabatic constant
gravity = 0.1 -- acceleration due to gravity
Lx = 1.0/6.0
Ly = 1.0
   
-- create app
eulerApp = Moments.App {

   tEnd = 8.5, -- end time
   nFrame = 5, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {Lx, Ly}, -- upper right corner
   cells = {50, 200}, -- number of cells
   cflFrac = 0.95, -- CFL fraction
   
   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { 
	 gasGamma = gasGamma,
      },

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local rhoTop = 2.0
	 local rhoBot = 1.0
	 local pr0 = 0.01
	 local pert = 0.01
	 local rho, pr
   
	 -- interface location
	 local yloc = 0.5 + pert*math.cos(2.0*math.pi*x/Lx)
	 if (y>yloc) then
	    -- in heavier fluid
	    rho = rhoTop
	    pr = rhoTop*gravity*(1-y)
	 else
	    -- in lighter fluid
	    rho = rhoBot
	    pr = rhoTop*gravity*(1-yloc) + rhoBot*gravity*(yloc-y)
	 end
   
	 return rho, 0.0, 0.0, 0.0, (pr+pr0)/(gasGamma-1)
      end,

      isAppliedAccelerationStatic = true,
      appliedAcceleration = function(t, xn)
	 return 0.0, -gravity, 0.0
      end,
      
      evolve = true, -- evolve species?

      bcx = { Moments.SpeciesBc.bcWall, Moments.SpeciesBc.bcWall },
      bcy = { Moments.SpeciesBc.bcWall, Moments.SpeciesBc.bcWall }
   }
}
-- run application
eulerApp:run()
