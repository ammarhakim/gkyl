-- Gkyl ------------------------------------------------------------------------
--
--

-- This test is just to check that the infrastructure supporting multiple
-- collision operators in the same species works.
-- The actual parameters may need to be adjusted later to make sure they
-- are consistent with some underlying assumptions in cross-collisions.

local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- Left/right state for shock.
local nl, ul, pl = 1.0, 0.0, 1.0
local nr, ur, pr = 0.125, 0.0, 0.1

local vthl = math.sqrt(pl/nl)
local vthr = math.sqrt(pr/nr)

local K = 0.005

local function maxwellian(n, u, vth, v)
   return n / math.sqrt(2*math.pi*vth^2) * 
      math.exp(-(v-u)^2 / (2*vth^2))
end

app = Plasma.App {
   logToFile = false,

   tEnd        = 0.015,             -- End time.
   nFrame      = 1,               -- Number of frames to write.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {1.0},            -- Configuration space upper right.
   cells       = {16},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 0.1,

   -- Decomposition for configuration space.
   decompCuts = {1},      -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {},     -- Periodic directions.

   neut1 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = {6.0},
      cells = {8},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

	 if math.abs(x) < 0.5 then
	    return maxwellian(nl, ul, vthl, v)
	 else
	    return maxwellian(nr, ur, vthr, v)
	 end
      end,
      bcx = { Plasma.Species.bcOpen,
	      Plasma.Species.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i" },
      -- Collisions.
      collBGK = Plasma.BgkCollisions {
      	 collideWith = { "neut1", "neut2", "neut3" },
      	 frequencies = { vthl/K, vthl/K, vthl/K }, 
      	 exactLagFixM012 = true,
      },
   },

   neut2 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = {6.0},
      cells = {8},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

	 if math.abs(x) < 0.5 then
	    return maxwellian(nl, ul, vthl, v)
	 else
	    return maxwellian(nr, ur, vthr, v)
	 end
      end,
      bcx = { Plasma.Species.bcOpen,
	      Plasma.Species.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i" },
     
      -- Collisions.
      collLBO = Plasma.LBOCollisions {
      	 collideWith = { "neut2", "neut3" },
      	 frequencies = { vthl/K/10, vthl/K/10 },
      },
      collBGK = Plasma.BgkCollisions {
      	 collideWith = { "neut1" },
      	 frequencies = { vthl/K }, 
      	 exactLagFixM012 = true,
      },
   },


   neut3 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = {6.0},
      cells = {8},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

	 if math.abs(x) < 0.5 then
	    return maxwellian(nl, ul, vthl, v)
	 else
	    return maxwellian(nr, ur, vthr, v)
	 end
      end,
      bcx = { Plasma.Species.bcOpen,
	      Plasma.Species.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i" },
     
      -- Collisions.
      collLBO = Plasma.LBOCollisions {
      	 collideWith = { "neut2", "neut3" },
      	 frequencies = { vthl/K/10, vthl/K/10 },
      },
      collBGK = Plasma.BgkCollisions {
      	 collideWith = { "neut1" },
      	 frequencies = { vthl/K }, 
      	 exactLagFixM012 = true,
      },
   },
}
-- Run application.
app:run()
