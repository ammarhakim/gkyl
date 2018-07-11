-- Gkyl ------------------------------------------------------------------------
local Hyper = require "App.HyperEqnOnCartGrid"
local Euler = require "Eq.Euler"
   
-- gas adiabatic index
gasGamma = 1.4

-- create app
eulerApp = Hyper.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {1024}, -- number of cells
   cfl = 0.9, -- CFL number
   limiter = "monotonized-centered", -- limiter
   equation = Euler { gasGamma = gasGamma }, -- equation to solve

   -- decomposition stuff
   decompCuts = {1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- initial condition
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
   
   -- boundary conditions
   periodicDirs = {}, -- periodic directions
   bcx = { Hyper.bcCopy, Hyper.bcCopy }, -- boundary conditions in X

   -- diagnostics
   diagnostics = {
      { name = "density", diagnostic = function (t, v) return v[1] end },
      { name = "totalEnergy", diagnostic = function (t, v) return v[5] end },
   },   
}
-- run appication
eulerApp:run()
