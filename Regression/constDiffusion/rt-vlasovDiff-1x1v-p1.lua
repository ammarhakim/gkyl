-- Gkyl -------------------------------------------------------------------
--
-- 1x1v simulation with diffusion along x or v using the Vlasov app,
-- but without evolving the collisionless terms.
--
---------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

local k_x = 0.3
local k_v = k_x

local nu_x = 0.1
local nu_v = nu_x

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 4./nu_x,          -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {-math.pi/k_x},   -- Configuration space lower left.
   upper       = { math.pi/k_x},   -- Configuration space upper right.
   cells       = {32},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".

   -- Boundary conditions for configuration space.
   periodicDirs = {1},   -- Periodic directions.
   
   -- Constant diffusion along x.
   fDiffx = Plasma.Species {
      charge = 0.,  mass = 1.,
      -- Velocity space grid.
      lower = {-math.pi/k_v},
      upper = { math.pi/k_v},
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         return (1+0.25*math.cos(k_x*x))
      end,
      evolve = true, -- Evolve species?
      evolveCollisionless = false, -- Don't evolve collisionless terms.
      -- Write out density, flow, total energy, and heat flux moments.
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      diff = Plasma.Diffusion {
         coefficient   = nu_x,
         -- Optional inputs:
         diffusiveDirs = {1},
         order = 2,   -- Diffusion order: 2, 4, or 6.
      },
   },

   -- Spatially varying diffusion along x.
   fVarDiffx = Plasma.Species {
      charge = 0.,  mass = 1.,
      -- Velocity space grid.
      lower = {-math.pi/k_v},
      upper = { math.pi/k_v},
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         return (1+0.25*math.cos(k_x*x))
      end,
      evolve = true, -- Evolve species?
      evolveCollisionless = false, -- Don't evolve collisionless terms.
      -- Write out density, flow, total energy, and heat flux moments.
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      diff = Plasma.Diffusion {
         coefficient   = {function(t,xn) return nu_x*(1+0.25*math.cos(2.*k_x*xn[1])) end},
         -- Optional inputs:
         diffusiveDirs = {1},
         order = 2,   -- Diffusion order: 2, 4, or 6.
      },
   },

}
-- Run application.
plasmaApp:run()
