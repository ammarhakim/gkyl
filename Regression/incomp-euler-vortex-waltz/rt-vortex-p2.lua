-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

local L = 10.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 100.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {0, 0}, -- configuration space lower left
   upper = {L, L}, -- configuration space upper right
   cells = {32,32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cfl = 0.2,

   -- decomposition for configuration space
   decompCuts = {1,1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   fluid = Plasma.IncompEuler.Species {
      charge = -1.0,
      -- initial conditions
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local x1, y1 = 3.5, 5.0
         local x2, y2 = 6.5, 5.0
         local r1 = (x-x1)^2 + (y-y1)^2
         local r2 = (x-x2)^2 + (y-y2)^2
         return math.exp(-r1/0.8^2) + math.exp(-r2/0.8^2) -- 4.0212385953656/L^2
      end,
      evolve = true, -- evolve species?
   },

   -- field solver
   field = Plasma.IncompEuler.Field {
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
