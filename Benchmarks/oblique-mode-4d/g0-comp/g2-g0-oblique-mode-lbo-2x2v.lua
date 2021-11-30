-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Constants
chargeElc = -1.0
massElc = 1.0

--R is ratio of vth to ud
ud = 0.3
R  = 0.333333333333333

-- Initial conditions
nElc10 = 0.5
nElc20 = 0.5
uxElc10 = 0.0
uyElc10 = ud
uzElc10 = 0.0
uxElc20 = 0.0
uyElc20 = -ud
uzElc20 = 0.0
TElc10 = massElc*(R*ud)^2
TElc20 = massElc*(R*ud)^2

k0 = 1.0
theta = 45.0/180.0*math.pi --0deg is pure Weibel, 90deg is pure Two Stream
kx = k0*math.cos(theta)
ky = k0*math.sin(theta)
perturb_n = 1e-8
alpha = 1.18281106421231 --ratio of E_y/E_x 

vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)

nuElc = 1.0e-4

Nx = 8 -- Number of configuration space point in x.
Ny = 8 -- Number of configuration space point in y.
Nvx = 16 -- Number of velocity space point in vx.
Nvy = 16 -- Number of velocity space point in vy.
Lx = 2*math.pi/kx -- domain size in x.
Ly = 2*math.pi/ky -- domain size in y.
Vmax = 0.9 -- Upper bound of velocity grid (3x initial drift velocity, 0.9 speed of light).

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = 15.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {0.0, 0.0}, -- configuration space lower left
   upper = {Lx, Ly}, -- configuration space upper right
   cells = {Nx, Ny}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {1,1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2}, -- periodic directions

   -- electrons 
   elc = Plasma.Species {
      charge = chargeElc, mass = massElc,
      -- velocity space grid
      lower = {-Vmax, -Vmax},
      upper = {Vmax, Vmax},
      cells = {Nvx, Nvy},
      -- initial conditions
      init = function (t, xn)
	 local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
	 local fv= maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
	    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)
	    return (1.0+perturb_n*math.cos(kx*x+ky*y))*fv
         end,
      evolve = true, -- evolve species?
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
      diagnosticMoments = { "M1i", "M2" }
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local E_x  = -perturb_n*math.sin(kx*x+ky*y)/(kx+ky*alpha)
	 local E_y  = alpha*E_x
         local B_z  = kx*E_y-ky*E_x
	 return E_x, E_y, 0.0, 0.0, 0.0, B_z
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
