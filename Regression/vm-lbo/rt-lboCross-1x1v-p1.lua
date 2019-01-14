-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

local n1, u1, p1 = 1.0, 0.0, 1.0
local n2, u2, p2 = 1.0, 1.0, 0.1

local vth1 = math.sqrt(p1/n1)
local vth2 = math.sqrt(p2/n2)

local K = 0.01

local function maxwellian(n, u, vth, v)
   return n / math.sqrt(2*math.pi*vth^2) * 
      math.exp(-(v-u)^2 / (2*vth^2))
end

print(" ")
print("Collisional period: ", K/vth1)
print("tEnd = ", 0.01)
print(" ")

vlasovApp = Plasma.App {
   logToFile = false,

   tEnd        = 0.01,              -- End time.
   nFrame      = 2,              -- Number of frames to write.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {1.0},            -- Configuration space upper right.
   cells       = {32},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- First neutral species.
   neut1 = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = { 6.0},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian(n1, u1, vth1, v)
      end,
      evolve = true,
      --evolveCollisions = false,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2" },

      -- Collisions.
      coll = Plasma.LBOCollisions {
	 collideWith  = { "neut1", "neut2", },
     	 frequencies  = { vth1/K, vth1/K, },
         -- Optional arguments:
         --crossOption = "Greene",    -- Or crossOption="HeavyIons".
         --betaGreene  = 1.0,
      },
   },

   -- Second neutral species.
   neut2 = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = { 6.0},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian(n2, u2, vth2, v)
      end,
      evolve = true,
      --evolveCollisions = false,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
	 collideWith = { "neut1", "neut2", },
     	 frequencies = { vth1/K, vth1/K, },
         -- Optional arguments:
         --crossOption = "Greene",    -- Or crossOption="HeavyIons".
         --betaGreene  = 1.0,
      },
   },
}
-- run application
vlasovApp:run()
