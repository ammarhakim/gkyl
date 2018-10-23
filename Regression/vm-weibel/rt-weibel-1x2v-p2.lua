-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell

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
k0 = 0.4
perturb = 1e-3
-- IC automatically calculated
vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)
vthIon0 = math.sqrt(TIon0/massIon)

knumber = 0.5 -- wave-number
elVTerm = 0.2 -- electron thermal velocity
vDrift = 1.0 -- drift velocity
perturbation = 1.0e-6 -- distribution function perturbation

-- Maxwellian in 1x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 50.0, -- end time
   nFrame = 1, -- number of output frames
   lower = { 0.0 }, -- configuration space lower left
   upper = { 2*math.pi/k0 }, -- configuration space upper right
   cells = {24}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Vlasov.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-1.0, -1.0},
      upper = {1.0, 1.0},
      cells = {12, 12},
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
   field = Vlasov.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x = xn[1]
	 local Bz = perturb*math.sin(k0*x)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, Bz
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
