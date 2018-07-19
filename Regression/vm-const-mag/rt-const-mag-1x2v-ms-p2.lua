-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- Constants
chargeElc = -1.0
chargeIon = 1.0
massElc = 1.0

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

vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)

knumber = 0.5 -- wave-number
elVTerm = 0.2 -- electron thermal velocity
vDrift = 1.0 -- drift velocity
B0 = 0.1 -- magnetic field

-- Maxwellian in 1x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = 2*math.pi/B0, -- end time
   nFrame = 1, -- number of output frames
   lower = { -1.0 }, -- configuration space lower left
   upper = { 1.0 }, -- configuration space upper right
   cells = {4}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Plasma.VlasovMaxwell.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-1.0, -1.0},
      upper = {1.0, 1.0},
      cells = {16, 16},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
	    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- field solver
   field = Plasma.VlasovMaxwell.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x = xn[1]
	 return 0.0, 0.0, 0.0, 0.0, 0.0, B0
      end,
      evolve = false, -- evolve field?
   },
}
-- run application
vlasovApp:run()
