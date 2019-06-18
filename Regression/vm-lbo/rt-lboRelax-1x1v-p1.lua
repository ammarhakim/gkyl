-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

polyOrder = 1
n0        = 1.0                             -- Density.
u0        = 0.0                             -- Flow speed.
vt        = 1.0/3.0                         -- Thermal speed..
nu        = 0.01                            -- Collision frequency.
-- The next three are for p1, v \in [-8vt,8vt], 2x32, rectangular IC.
nMr  = 1.0103629711+1.845e-12               -- Density of Maxwellian and rectangle. 
uMr  = 0.0                                  -- Flow speed of Maxwellian and rectangle. 
vtMr = math.sqrt(0.11235161663+1.62551e-12) -- Thermal speed of Maxwellian and rectangle.
-- Large bump on tail of Maxwellian:
ab   = math.sqrt(0.1)                       -- Amplitude of bump.
ub   = 4*math.sqrt(((3*vt/2)^2)/3)          -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0                                  -- Thermal speed of Maxwellian in bump.
-- The next three are for p1, v \in [-8vt,8vt], 2x32, bump in tail IC. 
nMb  = 1.1018707667+5.1e-13                 -- Density of Maxwellian and bump. 
uMb  = 0.4888097565+27.426e-12              -- Flow speed of Maxwellian and bump. 
vtMb = math.sqrt(0.39691067156+2.43e-13)    -- Thermal speed of Maxwellian and bump.

-- Top hat function without drift (u=0).
local function topHat(x, v, n, u, vth)
   local fOut   = 0.0
   local v0     = math.sqrt(3.0)*vth
   if math.abs(v) < v0 then
      fOut = n/(2.0*v0)
   else
      fOut = 0.0
   end
   return fOut
end

-- Maxwellian with a Maxwellian bump in the tail.
local function bumpMaxwell(x,vx,n,u,vth,bA,bU,bS,bVth)
   local Pi   = math.pi
   local vSq  = ((vx-u)/(math.sqrt(2.0)*vth))^2
   local vbSq = ((vx-u)/(math.sqrt(2.0)*bVth))^2

   return (n/math.sqrt(2.0*Pi*vth))*math.exp(-vSq)
         +(n/math.sqrt(2.0*Pi*bVth))*math.exp(-vbSq)*(bA^2)/((vx-bU)^2+bS^2)
end

plasmaApp = Plasma.App {
   logToFile = false,

   tEnd        = 100,           -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower coordinate.
   upper       = {1.0},         -- Configuration space upper coordinate.
   cells       = {2},           -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,     -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space.
   decompCuts = {1},            -- Cuts in each configuration direction.
   useShared  = false,          -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},          -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   square = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower      = {-8.0*vt},
      upper      = {8.0*vt},
      cells      = {48},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

         return topHat(x, v, n0, u0, vt)
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'square'},
         frequencies = {nu},
      },
   },

   -- -- Maxwellian for comparison with rectangular IC.
   -- maxwellSquare = Plasma.Species {
   --    charge = 0.0, mass = 1.0,
   --    -- Velocity space grid.
   --    lower      = {-8.0*vt},
   --    upper      = {8.0*vt},
   --    cells      = {48},
   --    -- Initial conditions.
   --    init = Plasma.MaxwellianProjection {
   --       density         = nMr,
   --       driftSpeed      = {uMr},
   --       temperature     = vtMr^2,
   --       exactScaleM0    = false,
   --       exactLagFixM012 = true,
   --    },
   --    -- Evolve species?
   --    evolve = false,
   --    -- Diagnostic moments.
   --    diagnosticMoments = { "M0", "M1i", "M2" },
   -- },

   -- Neutral species with a bump in the tail.
   bump = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower      = {-8.0*vt},
      upper      = {8.0*vt},
      cells      = {48},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

         return bumpMaxwell(x,v,n0,u0,vt,ab,ub,sb,vtb)
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'bump'},
         frequencies = {nu},
      },
   },

   -- -- Maxwellian for comparison with bump in tail IC.
   -- maxwellBump = Plasma.Species {
   --    charge = 0.0, mass = 1.0,
   --    -- Velocity space grid.
   --    lower      = {-8.0*vt},
   --    upper      = {8.0*vt},
   --    cells      = {48},
   --    -- Initial conditions.
   --    init = Plasma.MaxwellianProjection {
   --       density         = nMb,
   --       driftSpeed      = {uMb},
   --       temperature     = vtMb^2,
   --       exactScaleM0    = false,
   --       exactLagFixM012 = true,
   --    },
   --    -- Evolve species?
   --    evolve = false,
   --    -- Diagnostic moments.
   --    diagnosticMoments = { "M0", "M1i", "M2" },
   -- },

}
-- run application
plasmaApp:run()
