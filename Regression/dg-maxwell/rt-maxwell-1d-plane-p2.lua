-- Gkyl ------------------------------------------------------------------------
--
-- 

local Vlasov    = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Constants = require "Lib.Constants"

-- Physical parameters.
epsilon0   = Constants.EPSILON0
mu0        = Constants.MU0
lightSpeed = Constants.SPEED_OF_LIGHT

L       = 1.0
kwave   = 2
freq    = 2*math.pi/L*math.sqrt(kwave^2)*lightSpeed
tperiod = 2*math.pi/freq

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd        = tperiod,       -- End time.
   nFrame      = 1,             -- Number of output frames.
   lower       = {-L/2.},       -- Configuration space lower left.
   upper       = { L/2.},       -- Configuration space upper right.
   cells       = {32},          -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Field solver.
   field = Vlasov.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      elcErrorSpeedFactor = 0,
      mgnErrorSpeedFactor = 0,
      
      init = function (t, xn)
	 local x, y  = xn[1], 0.0
	 local cos   = math.cos
	 local pi    = math.pi
	 local c     = lightSpeed
	 local phi   = 2*pi/L*(kwave*x)
	 local knorm = math.sqrt(kwave^2)
	 local kxn   = kwave/knorm
         -- n = (0, 1, 1), n_hat = 1/math.sqrt(2)
	 local E0 = 1.0/math.sqrt(2.0)
	 local Ex = 0.0
         local Ey = E0*cos(phi)
	 local Ez = E0*cos(phi)
	 local Bx = 0.0
	 local By = -E0/c*cos(phi)*kxn
	 local Bz = E0/c*cos(phi)*kxn
	 return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- Evolve field?
   },   
   
}
-- Run application.
vlasovApp:run()
