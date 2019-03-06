-- Plasma ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").IncompEuler

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

   tEnd        = 1.0,              -- End time.
   nFrame      = 20,               -- Number of output frames.
   lower       = {0, 0},           -- Configuration space lower left.
   upper       = {1.0, 1.0},       -- Configuration space upper right.
   cells       = {16, 16},         -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   
   -- Decomposition for configuration space.
   decompCuts = {1,1},    -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2}, -- Periodic directions

   -- Fluid species.
   fluid = Plasma.Species {
      charge = 1.0,
      -- Initial conditions.
      init            = squareHat,
      evolve          = true, -- Evolve species?
      applyPositivity = true,
      diff = Plasma.Diffusion {
         coefficient = {0.1, 0.1},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = false, -- Evolve field?
      -- u = {dphi/dy, -dphi/dx}
      initPhiFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         return -uy*x + ux*y    -- Moving in \xhat+\yhat direction. 
--         return  ux*y           -- Moving in \xhat direction. 
--         return -uy*x           -- Moving in \yhat direction. 
      end, 
   },
}
-- Run application.
plasmaApp:run()
