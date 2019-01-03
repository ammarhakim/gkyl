-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

local ux = 1
local uy = 1

local function gaussian(t, xn)
   local r2 = (xn[1]-0.5)^2 + (xn[2]-0.5)^2
   return math.exp(-50*r2)
end
local function cylinder(t, xn)
   local r2 = (xn[1]-0.5)^2 + (xn[2]-0.5)^2
   if r2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
local function step(t, xn)
   local r2 = (xn[1]-0.5)^2
   if r2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
local function squareHat(t, xn)
   local rx2, ry2 = (xn[1]-0.5)^2, (xn[2]-0.5)^2
   if rx2 < 0.25^2 and ry2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
local function expTent(t, xn)
   local r = math.sqrt((xn[1]-0.5)^2 + (xn[2]-0.5)^2)
   return math.exp(-10*r)
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = .10, -- end time
   nFrame = 10, -- number of output frames
   lower = {0, 0}, -- configuration space lower left
   upper = {1.0, 1.0}, -- configuration space upper right
   cells = {16, 16}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   
   -- decomposition for configuration space
   decompCuts = {1,1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   fluid = Plasma.IncompEuler.Species {
      charge = 1.0,
      -- initial conditions
      init = squareHat,
      evolve = true, -- evolve species?
      applyPositivity = true,
   },

   -- field solver
   field = Plasma.IncompEuler.Field {
      evolve = false, -- evolve field?
      -- u = {dphi/dy, -dphi/dx}
      initPhiFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         return -uy*x + ux*y 
      end, 
   },
}
-- run application
plasmaApp:run()
