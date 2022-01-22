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
-- Large bump on tail of Maxwellian:
vt   = math.sqrt(1.0/24.0)                  -- Thermal speed of Maxwellian in bump.
ab   = 4*math.sqrt(0.1)                     -- Amplitude of bump.
ub   = {2*math.sqrt(1/3),0.0}               -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0/math.sqrt(2.0)                   -- Thermal speed of Maxwellian in bump.

vLower0 = {-2.0,-2.0}                     -- Min velocity in grid.
vUpper0 = { 2.0, 2.0}                     -- Max velocity in grid.
Nx      = {2}                             -- Number of cells in configuration space.
Nv0     = {16,16}                         -- Number of cells in velocity space.

extraCells = {6,6}
dv     = {(vUpper0[1]-vLower0[1])/Nv0[1], (vUpper0[2]-vLower0[2])/Nv0[2]}
vLower = {vLower0[1]-(extraCells[1]/2)*dv[1], vLower0[2]-(extraCells[2]/2)*dv[2]}
vUpper = {vUpper0[1]+(extraCells[1]/2)*dv[1], vUpper0[2]+(extraCells[2]/2)*dv[2]}
Nv     = {Nv0[1]+extraCells[1], Nv0[2]+extraCells[2]}

finalTime = 25.0 
numFrames = 1

maskOut = function(t,xn)
   local x, vx, vy = xn[1], xn[2], xn[3]
   if vx < vLower0[1] or vx > vUpper0[1]
      or vy < vLower0[2] or vy > vUpper0[2] then
      return -1
   else
      return 1
   end
end

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

   tEnd        = finalTime,     -- End time.
   nFrame      = numFrames,     -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower coordinate.
   upper       = {1.0},         -- Configuration space upper coordinate.
   cells       = Nx,       -- Configuration space cells.
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
      lower = vLower,
      upper = vUpper,
      cells = Nv,
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         return topHat(x, vx, vy, n0, u0[1], u0[2], vt)
      end,
      mask = maskOut,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
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
      lower = vLower,
      upper = vUpper,
      cells = Nv,
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         return bumpMaxwell(x,vx,vy,n0,u0[1],u0[2],vt,ab,ub[1],ub[2],sb,vtb)
      end,
      mask = maskOut,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'bump'},
         frequencies = {nu},
      },
   },

}
-- Run application.
plasmaApp:run()
