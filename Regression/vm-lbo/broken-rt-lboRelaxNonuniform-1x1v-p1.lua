-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

polyOrder = 1
n0        = 1.0                             -- Density.
u0        = 0.0                             -- Flow speed.
vt        = 1.0/3.0                         -- Thermal speed..
nu        = 0.01                            -- Collision frequency.
-- Large bump on tail of Maxwellian:
ab   = math.sqrt(0.1)                       -- Amplitude of bump.
ub   = 4*math.sqrt(((3*vt/2)^2)/3)          -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0                                  -- Thermal speed of Maxwellian in bump.

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

   tEnd        = 20,            -- End time.
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
      lower = {-1.0},
      upper = { 1.0},
      coordinateMap = {
        function(z) if z<0. then return -8.0*vt*math.abs(z)^2 else return 8.0*vt*z^2 end end
      },
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

         return topHat(x, v, n0, u0, vt)
      end,
      -- Evolve species?
      evolve = true,
      diagnostics = { "M0", "M1i", "M2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'square'},
         frequencies = {nu},
      },
   },

   -- Neutral species with a bump in the tail.
   bump = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-1.0},
      upper = { 1.0},
      coordinateMap = {
        function(z) if z<0. then return -8.0*vt*math.abs(z)^2 else return 8.0*vt*z^2 end end
      },
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]

         return bumpMaxwell(x,v,n0,u0,vt,ab,ub,sb,vtb)
      end,
      -- Evolve species?
      evolve = true,
      diagnostics = { "M0", "M1i", "M2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'bump'},
         frequencies = {nu},
      },
   },

}
-- Run application.
plasmaApp:run()
