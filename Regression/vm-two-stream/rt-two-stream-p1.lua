-- Gkyl ------------------------------------------------------------------------
local Vlasov = require "App.VlasovOnCartGrid"

knumber = 0.5 -- wave-number
elVTerm = 0.2 -- electron thermal velocity
vDrift = 1.0 -- drift velocity
perturbation = 1.0e-6 -- distribution function perturbation

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 40.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
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
      lower = {-6.0},
      upper = {6.0},
      cells = {32},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local alpha = perturbation
	 local k = knumber
	 local vt = elVTerm
	 
	 local fv = 1/math.sqrt(8*math.pi*vt^2)*(math.exp(-(v-vDrift )^2/(2*vt^2))+math.exp(-(v+vDrift)^2/(2*vt^2)))
	 return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- field solver
   field = Vlasov.EmField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local alpha = perturbation
	 local k = knumber
	 return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
