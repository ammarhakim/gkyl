-- Gkyl ------------------------------------------------------------------------
--
-- 1x1v p=2 ion acoustic instability simulation.
--
-------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Electron parameters.
vDriftElc = {0.01,0}   -- Modified from 0.159.
vtElc     = 0.02
-- Ion parameters.
vDriftIon = {0.0,0.0}
vtIon     = 0.000566      -- Modified from 0.001 (use Te = 50 Ti).
-- Mass ratio.
massRatio = 25.0  -- Modified from 25. 
tempRatio = 50.0
nuee = 0.001
nuii = nuee/math.sqrt(massRatio)*math.pow(tempRatio,1.5)

knumberx     = 10.0   -- Wave-number.
knumbery     = 0.0
perturbation = 1.0e-4 -- Distribution function perturbation.

local function maxwellian2v(v, vDrift, vt)
    return 1/(2*math.pi*vt^2)*math.exp(-((v[1]-vDrift[1])^2+(v[2]-vDrift[2])^2)/(2*vt^2))
end

plasmaApp = Plasma.App {
   logToFile = true,
 
   tEnd        = 2.0,             -- End time.
   nFrame      = 1,               -- Number of output frames.
   nDistFuncFrame = 1,            -- Number of distribution function output frames 

   lower       = {0.0,0.0},        -- Configuration space lower left.
   upper       = {1.0,1.0},        -- Configuration space upper right.
   cells       = {16,16},          -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 0.9,
 
   -- Decomposition for configuration space.
   decompCuts = {1, 1},      -- Cuts in each configuration direction.
   useShared  = false,       -- If to use shared memory.
 
   -- Boundary conditions for configuration space.
   periodicDirs = {1,2}, -- Periodic directions.
   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.01,
 
   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-3.0*vtElc,-3.0*vtElc},
      upper = { 3.0*vtElc, 3.0*vtElc},
      cells = {24,24},
      decompCuts = {1,1},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
         local fv = maxwellian2v({vx, vy}, vDriftElc, vtElc)
         return fv 
      end,
      evolve = true, -- Evolve species?

      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuee},
      },
 
      diagnostics = { "M0", "M1i", "intM0", "intM1i" },
   },
 
   -- Ions.
   ion = Plasma.Species {
      charge = 1.0, mass = massRatio,
      -- Velocity space grid.
      lower = {-6.0*vtIon,-6.0*vtIon},
      upper = { 6.0*vtIon, 6.0*vtIon},
      cells = {24,24},
      decompCuts = {1,1},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
         local fv = maxwellian2v({vx, vy}, vDriftIon, vtIon)
         return fv*(1+perturbation*math.cos(2*math.pi*knumberx*x+ 2*math.pi*knumbery*y))
      end,
      evolve = true,    -- Evolve species?

      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuii},
      },
 
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0,
      evolve   = true, -- Evolve field?
      hasMagneticField = false,
   },
 
}
-- Run application.
plasmaApp:run()
 
