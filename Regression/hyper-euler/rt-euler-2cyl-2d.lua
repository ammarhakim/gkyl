-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 1.4 -- gas adiabatic constant

function circle(x, y, x0, y0, rad)
   return (x-x0)^2+(y-y0)^2<rad^2
end

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.16, -- end time
   nFrame = 5, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {1.0, 1.0}, -- upper right corner
   cells = {300, 300}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux="lax" },
      forceInv = false,

      -- initial conditions
      init = function (t, xn)
         local x = xn[1]
         -- Mach 2 shock
         local rho, u, pr =   3.85714285*1.4,  2.22222222, 10.3333333*1.0
         if (x>0.2) then
            rho, u, pr = 1.4, 0.0, 1.0
         end
         return rho, rho*u, 0, 0, pr/(gasGamma-1) + 0.5*rho*u*u
      end,
      evolve = true, -- evolve species?

      -- outer boundary conditions
      bcx = { Moments.Species.bcCopy,
              Moments.Species.bcCopy },
      bcy = { Moments.Species.bcCopy, Moments.Species.bcCopy },

      -- has interior (embedded) boundary conditions?
      hasSsBnd = true,
      -- mask defining the embedded boundary
      inOutFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         local xc, yc = 0.0, 0.0
         local rad = 0.5
         return  (circle(x,y, 0.4, 0.25, 0.15) or circle(x,y, 0.5, 0.75, 0.15)) and -1.0 or 1.0
      end,
      -- boundary conditions to be applied on the embedded boundary
      ssBc = { Moments.Species.bcWall },
   },   
}

-- run application
eulerApp:run()
