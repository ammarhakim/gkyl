-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- Constants
chargeElc = -1.0
chargeIon = 1.0
massElc = 1.0
massIon = 1836.0

-- Initial conditions
nElc10 = 0.5
nElc20 = 0.5
uxElc10 = 0.0
uyElc10 = 0.3
uzElc10 = 0.0
uxElc20 = 0.0
uyElc20 = -0.3
uzElc20 = 0.0
TElc10 = 0.01
TElc20 = 0.01
nIon0 = 1.0
uxIon0 = 0.0
uyIon0 = 0.0
uzIon0 = 0.0
TIon0 = 0.01
k0 = 1.0
theta = 45.0/180.0*math.pi --0deg is pure Weibel, 90deg is pure Two Stream
kx = k0*math.cos(theta)
ky = k0*math.sin(theta)
perturb_n = 1e-4
alpha = 0.23669277173847925 --ratio E_y/E_x (angle and k0 specific!)

vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)
vthIon0 = math.sqrt(TIon0/massIon)

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = 30.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {0, 0}, -- configuration space lower left
   upper = {2*math.pi/kx, 2*math.pi/ky }, -- configuration space upper right
   cells = {8, 8}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1,1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2}, -- periodic directions

   -- electrons 
   elc = Plasma.VlasovMaxwell.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-1.0, -1.0},
      upper = {1.0, 1.0},
      cells = {8, 8},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
	 local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
	 local fv= maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
	    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)
	 return fv*(1.0+perturb_n*math.cos(kx*x+ky*y))
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- field solver
   field = Plasma.VlasovMaxwell.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local E_x  = -perturb_n*math.sin(kx*x+ky*y)/kx
	 local E_y  = -alpha*perturb_n*math.sin(kx*x+ky*y)/ky
	 local Bz   = kx*E_y-ky*E_x
	 return E_x, E_y, 0.0, 0.0, 0.0, Bz
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
