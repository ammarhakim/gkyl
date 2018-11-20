-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 4*math.pi, -- end time
   nFrame = 1, -- number of output frames
   lower = {0, 0}, -- configuration space lower left
   upper = {1.0, 1.0}, -- configuration space upper right
   cells = {32, 32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   
   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   fluid = Plasma.IncompEuler.Species {
      charge = 1.0,
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local x0, y0, r0 = 0.25, 0.5, 0.15
	 local r = math.min(math.sqrt((x-x0)^2+(y-y0)^2), r0)/r0
	 return 0.25*(1+math.cos(math.pi*r))
      end,
      evolve = true, -- evolve species?
      applyPositivity = true,
   },

   -- field solver
   field = Plasma.IncompEuler.Field {
      evolve = false, -- evolve field?
      -- u = {dphi/dy, -dphi/dx}
      initPhiFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         return -0.5*(y^2-y+x^2-x)
      end, 
   },
}
-- run application
plasmaApp:run()
