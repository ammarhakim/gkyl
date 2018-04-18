-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

TiTe = 1.0
knumber = 0.5
XL, XU = -math.pi/knumber, math.pi/knumber
Ti0 = 1.0
Te0 = Ti0/TiTe
ni0, ne0 = 1.0, 1.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 15, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- gyrokinetic ions
   ion = Plasma.GkSpecies {
      charge = 1.0,
      mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {128},
      decompCuts = {1},
      -- initial conditions
      -- specify background so that we can plot perturbed distribution and moments
      initBackground = function (t, xn) 
	 local x, v = xn[1], xn[2]
         local k = knumber
         local vth2 = Ti0
         local alpha = 0.01
	 return ni0/math.sqrt(2*math.pi*vth2)*math.exp(-v^2/(2*vth2))
      end,
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local k = knumber
         local vth2 = Ti0
         local alpha = 0.01
	 return ni0/math.sqrt(2*math.pi*vth2)*math.exp(-v^2/(2*vth2))*(1.0 + alpha * math.cos(k*x))
      end,
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkDens"},
   },

   adiabaticElectron = Plasma.GkSpecies {
      charge = -1.0,
      mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {128},
      decompCuts = {1},
      -- initial conditions.. use ion background so that background is exactly neutral
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local k = knumber
         local vth2 = Ti0
	 return ni0/math.sqrt(2*math.pi*vth2)*math.exp(-v^2/(2*vth2))
      end,
      evolve = false, -- evolve species?
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve field?
      adiabatic = {response = "electron", charge = -1.0, dens = ne0, temp = Te0},
      discontinuous = false,
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         return 1.0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
