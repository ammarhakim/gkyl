-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Constants = require "Lib.Constants"

X = 80.0 -- [m]
Y = 40.0 -- [m]
NX = 80
NY = 40
   
-- create app
maxwellApp = Moments.App {
   logToFile = true,

   tEnd = 150e-9, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {X, Y}, -- upper right corner
   cells = {NX, NY}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {}, -- periodic directions

   -- electrons
   field = Moments.Field {
      epsilon0 = Constants.EPSILON0, mu0 = Constants.MU0,
      mgnErrorSpeedFactor = 1.0,

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local m, n = 8, 5
	 local a = m*math.pi/X
	 local b = n*math.pi/Y
	 local Ez = math.sin(a*x)*math.sin(b*y)
	 return 0.0, 0.0, Ez, 0.0, 0.0, 0.0
      end,

      evolve = true, -- evolve species?

      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect},
      bcy = { Moments.Field.bcReflect, Moments.Field.bcReflect},
   },   
}
-- run application
maxwellApp:run()
