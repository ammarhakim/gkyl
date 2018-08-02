-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Constants = require "Lib.Constants"

L = 1.0
X = L -- [m]
Y = L -- [m]
NX = 40
NY = 40
kwave = 2
lwave = 2
freq = 2*math.pi/L*math.sqrt(kwave^2+lwave^2)*Constants.SPEED_OF_LIGHT

-- create app
maxwellApp = Moments.App {
   logToFile = true,

   tEnd = 2*math.pi/freq, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {X, Y}, -- upper right corner
   cells = {NX, NY}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   field = Moments.Field {
      epsilon0 = Constants.EPSILON0, mu0 = Constants.MU0,
      mgnErrorSpeedFactor = 1.0,

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local c = Constants.SPEED_OF_LIGHT
	 local phi = 2*math.pi/L*(kwave*x + lwave*y)
	 local E0 = 1.0
	 local Ex, Ey = 0.0, 0.0
	 local Ez = E0*math.cos(phi)
	 local Bx = E0/c*math.cos(phi)/math.sqrt(2)
	 local By = -E0/c*math.cos(phi)/math.sqrt(2)
	 local Bz = 0.0
	 return Ex, Ey, Ez, Bx, By, Bz
      end,

      evolve = true, -- evolve species?
   },   
}
-- run application
maxwellApp:run()
