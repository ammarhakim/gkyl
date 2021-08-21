-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

polyOrder = 2
nu        = 0.01                            -- Collision frequency.
n0        = 1.0                             -- Density.
u0        = {0.0,0.0}                       -- Flow speed.
vMin      = -2.0                            -- Min velocity in grid.
vMax      =  2.0                            -- Max velocity in grid.
Nx        = {2}                             -- Number of cells in configuration space.
Nv        = {16,16}                         -- Number of cells in velocity space.
-- Large bump on tail of Maxwellian:
vt   = math.sqrt(1.0/24.0)                  -- Thermal speed of Maxwellian in bump.
ab   = 4*math.sqrt(0.1)                     -- Amplitude of bump.
ub   = {2*math.sqrt(1/3),0.0}               -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0/math.sqrt(2.0)                   -- Thermal speed of Maxwellian in bump.

-- Top hat function without drift (u=0).
local function topHat(x, vx, vy, n, ux, uy, vth)
   local fOut   = 0.0
   local v0     = math.sqrt(3.0*(vth^2)/2.0)
   if (math.abs(vx) <= v0) and (math.abs(vy) <= v0) then
      fOut = n/((2.0*v0)^2)
   else
      fOut = 0.0
   end
   return fOut
end

-- Maxwellian with a Maxwellian bump in the tail.
local function bumpMaxwell(x,vx,vy,n,ux,uy,vth,bA,bUx,bUy,bS,bVth)
   local Pi   = math.pi
   local vSq  = ((vx-ux)^2+(vy-uy)^2)/((math.sqrt(2.0)*vth)^2)
   local vbSq = ((vx-ux)^2+(vy-uy)^2)/((math.sqrt(2.0)*bVth)^2)

   return (n/(2.0*Pi*(vth^2)))*math.exp(-vSq)
         +(n/(2.0*Pi*(bVth^2)))*math.exp(-vbSq)*(bA^2)/((vx-bUx)^2+(vy-bUy)^2+bS^2)
end

plasmaApp = Plasma.App {
   logToFile = false,

   tEnd        = 50,           -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower coordinate.
   upper       = {1.0},         -- Configuration space upper coordinate.
   cells       = {Nx[1]},       -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,     -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".
   -- unnecessary, but linear trigger bug causes multiple write outs at end 
   -- without small adjustment of time step
   cflFrac     = 0.95,

   -- Decomposition for configuration space.
   decompCuts = {1},            -- Cuts in each configuration direction.
   useShared  = false,          -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},          -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   square = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {vMin,vMin},
      upper = {vMax,vMax},
      cells = Nv,
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         return topHat(x, vx, vy, n0, u0[1], u0[2], vt)
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
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
      lower = {vMin,vMin},
      upper = {vMax,vMax},
      cells = Nv,
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         return bumpMaxwell(x,vx,vy,n0,u0[1],u0[2],vt,ab,ub[1],ub[2],sb,vtb)
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
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
