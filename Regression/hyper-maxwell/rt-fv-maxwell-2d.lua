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

   tEnd = 200e-9, -- end time
   nFrame = 5, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {X, Y}, -- upper right corner
   cells = {NX, NY}, -- number of cells
   timeStepper = "fvDimSplit",
   cflFrac = 0.5,
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {2}, -- periodic directions

   -- electrons
   field = Moments.Field {
      epsilon0 = Constants.EPSILON0, mu0 = Constants.MU0,
      mgnErrorSpeedFactor = 0.0,
      -- limiter = "no-limiter",
      limiter = "zero",

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local m, n = 8, 5
	 local a = m*math.pi/X
	 local Ez = math.sin(a*x)
	 local By = math.cos(a*x)/Constants.SPEED_OF_LIGHT
	 return 0.0, 0.0, Ez, 0.0, By, 0.0
      end,

      evolve = true, -- evolve species?

      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect},
   },   
}
-- run application
maxwellApp:run()
