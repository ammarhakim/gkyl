-- Gkyl ------------------------------------------------------------------------
--
-- Perform the test in vanLeer's 2005 recovery paper of a diffusion equation
-- with a source.
--
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").IncompEuler()

local diffCoefficient = 0.01

-- Options to benchmark against Van Leer (2005) using projection.
local function vanLeerIC(t, xn)
   local x = xn[1]
   local k = 2.0*math.pi
   return math.cos(k*x)
end
local function vanLeerSource(t, xn)
   local x = xn[1]
   local k = 2.0*math.pi
   return 4.0*(math.pi^2)*diffCoefficient*math.sin(k*x)
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 100.0,            -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {1.0},            -- Configuration space upper right.
   cells       = {16},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   
   -- Decomposition for configuration space.
   decompCuts = {1},      -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Fluid species.
   fluid = Plasma.Species {
      charge = 1.0,
      -- Projected initial conditions (with quadrature).
      init   = vanLeerIC,
      source = Plasma.Source{
         profile        = vanLeerSource,
         timeDependence = function (t) return 1.0 end,
      },
      evolve              = true, -- Evolve species?
      evolveCollisionless = false,
      diff = Plasma.Diffusion {
         coefficient = diffCoefficient,
      },
      -- BCs for the vanLeer 2005 test.
      bcx = {Plasma.BasicBC{kind="dirichlet", value=1.0},
             Plasma.BasicBC{kind="neumann", value=2.0*math.pi-1.0}},
   },

}
-- Run application.
plasmaApp:run()
