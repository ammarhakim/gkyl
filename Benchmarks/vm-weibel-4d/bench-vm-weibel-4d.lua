-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local prng = require "sci.prng"
local rng = prng.mrg32k3a()

-- Constants
chargeElc = -1.0
chargeIon = 1.0
massElc = 1.0
massIon = 1836.0

--R is ratio of vth to ud
ud = 0.1
R  = 0.1

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
nIon0 = 1.0
uxIon0 = 0.0
uyIon0 = 0.0
uzIon0 = 0.0
TIon0 = 0.01

k0_TS     = 6.135907273413176
theta = 90.0/180.0*math.pi --0deg is pure Weibel, 90deg is pure Two Stream
kx_TS = k0_TS*math.cos(theta)
ky_TS = k0_TS*math.sin(theta)

k0_Weibel = 2.31012970008316
theta = 0.0/180.0*math.pi --0deg is pure Weibel, 90deg is pure Two Stream
kx_Weibel = k0_Weibel*math.cos(theta)
ky_Weibel = k0_Weibel*math.sin(theta)
kx    = k0_Weibel
ky    = k0_TS/3.0	

perturb_n = 1e-4*ky_TS
perturb_B = perturb_n/ky_TS--1e-3
N=16
P={}
for i=0,N,1 do
   P[i]={}
   for j=0,N,1 do
      P[i][j]={}
      for k=1,6,1 do         
	 P[i][j][k]=rng:sample()
      end
   end
end

vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)
vthIon0  = math.sqrt(TIon0/massIon)

nuElc = 1.0e-5

Lx = 2*math.pi/kx
Ly = 2*math.pi/ky
Nx = 32
Ny = 32

vLimElc = 3*ud
NvElc = 32

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = 50.0, -- end time
   nFrame = 1, -- number of output frames
   lower = { 0, 0 }, -- configuration space lower left
   upper = { Lx, Ly }, -- configuration space upper right

   cells = {Nx,Ny}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3s4", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1,1}, -- cuts in each configuration direction
   useShared = true, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2}, -- periodic directions

   --restartFrameEvery = 0.02,
   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,

   -- electrons
   elc = Plasma.Species {
      charge = chargeElc, mass = massElc,
      -- velocity space grid
      lower = {-vLimElc, -vLimElc},
      upper = {vLimElc, vLimElc},
      cells = {NvElc, NvElc},
      -- initial conditions
      init = function (t, xn)
         local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
         local fv = maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
            maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)
	 return fv
      end,
      evolve = true, -- evolve species?
      diagnosticMoments = {"M0","M1i","M2ij","M3i"},
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x, y, E_x, E_y, B_z = xn[1], xn[2], 0.0, 0.0, 0.0
	 for i=0,N,1 do
	    for j=0,N,1 do
               if i~=0 or j~=0 then          
		  E_x  = E_x + perturb_n/N/N*P[i][j][1]*math.sin(i*kx*x+j*ky*y+2*math.pi*P[i][j][2])
		  E_y  = E_y + perturb_n/N/N*P[i][j][3]*math.sin(i*kx*x+j*ky*y+2*math.pi*P[i][j][4])
		  B_z  = B_z + perturb_n/N/N*P[i][j][5]*math.sin(i*kx*x+j*ky*y+2*math.pi*P[i][j][6])
               end
	    end
         end
         return E_x, E_y, 0.0, 0.0, 0.0, B_z
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
