-- Gkyl ------------------------------------------------------------------------
--
-- 

local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell
local Constants = require "Lib.Constants"

-- Physical parameters
epsilon0 = Constants.EPSILON0
mu0 = Constants.MU0
lightSpeed = Constants.SPEED_OF_LIGHT

L = 1.0
-- require k dot n = 0 so we choose n = (-1, 1, 1)
kwave = 3
lwave = 2
mwave = 1
freq = 2*math.pi/L*math.sqrt(kwave^2+lwave^2+mwave^2)*lightSpeed
tperiod = 2*math.pi/freq

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = tperiod, -- end time
   nFrame = 1, -- number of output frames
   lower = {0.0, 0.0, 0.0}, -- configuration space lower left
   upper = {L, L, L}, -- configuration space upper right
   cells = {16, 16, 16}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {1, 1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2, 3}, -- periodic directions

   -- field solver
   field = Vlasov.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      elcErrorSpeedFactor = 0,
      mgnErrorSpeedFactor = 0,
      
      init = function (t, xn)
	 local x, y, z = xn[1], xn[2], xn[3]
	 local cos = math.cos
	 local pi = math.pi
	 local c = lightSpeed
         local phi = 2*pi/L*(kwave*x+lwave*y+mwave*z)
	 local knorm = math.sqrt(kwave^2+lwave^2+mwave^2)
	 local kxn, kyn, kzn = kwave/knorm, lwave/knorm, mwave/knorm
         -- n_hat = 1/math.sqrt(3.0) if n = (-1,1,1)
	 local E0 = 1.0/math.sqrt(3.0)
         local Ex = -E0*cos(phi)
         local Ey = E0*cos(phi)
         local Ez = E0*cos(phi)
         -- some k x n have sign flips due to n = (-1,1,1)
         local Bx = E0/c*cos(phi)*2*pi/L*(kyn - kzn)
         local By = E0/c*cos(phi)*2*pi/L*(-kzn - kxn)
         local Bz = E0/c*cos(phi)*2*pi/L*(kxn + kyn)
	 return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- evolve field?
   },   
   
}
-- run application
vlasovApp:run()
