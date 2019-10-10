-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

polyOrder = 2
nu        = 0.01                            -- Collision frequency.
n0        = 1.0                             -- Density.
u0        = 0.0                             -- Flow speed.
B0        = 1.0                             -- Background magnetic field
mass      = 1.0                             -- Define mass globally
vMin      = -2.0                            -- Min velocity in grid.
vMax      =  2.0                            -- Max velocity in grid.
muMin     = 0.0 
muMax     = mass*(vMax^2)/(2*B0)
Nx        = {2,2,2}                             -- Number of cells in configuration space.
Nv        = {16,16}                         -- Number of cells in velocity space.
-- Large bump on tail of Maxwellian:
vt   = math.sqrt(1.0/12.0)                  -- Thermal speed of Maxwellian in bump.
ab   = math.sqrt(0.1)                       -- Amplitude of bump.
ub   = 4*math.sqrt(1/12.0)                    -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0/math.sqrt(2.0)                   -- Thermal speed of Maxwellian in bump.

-- Top hat function without drift (u=0).
local function topHat(x, y, z, vpar, mu, n, u, vth)
   local Pi     = math.pi
   local fOut   = 0.0
   local v0     = math.sqrt(3.0*(vth^2)/2.0)
   if (math.abs(vpar) <= v0) and (mu <= v0^2*mass/(2*B0)) then
      fOut = (1.0/(2.0*Pi)/v0^3)
   else
      fOut = 0.0
   end
   return fOut
end

-- Maxwellian with a Maxwellian bump in the tail.
local function bumpMaxwell(x,y,z,vpar,mu,n,u,vth,bA,bU,bS,bVth)
   local Pi   = math.pi
   local vSq  = ((vpar-u)^2+2*B0*math.abs(mu)/mass)/((math.sqrt(2.0)*vth)^2)
   local vbSq = ((vpar-u)^2+2*B0*math.abs(mu)/mass)/((math.sqrt(2.0)*bVth)^2)

   return (n/math.sqrt(2.0*Pi*(vth^2))^(3/2))*math.exp(-vSq)
      +(n/math.sqrt(2.0*Pi*(bVth^2))^(3/2))*math.exp(-vbSq)*(bA^2)/((vpar-bU)^2+bS^2)
end

plasmaApp = Plasma.App {
   logToFile = false,
   
   tEnd        = 0.01,           -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0,0,0},       -- Configuration space lower coordinate.
   upper       = {1.0,1.0,1.0}, -- Configuration space upper coordinate.
   cells       = Nx,            -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,     -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 0.01,

   
   -- Decomposition for configuration space.
   decompCuts = {1,1,1},            -- Cuts in each configuration direction.
   useShared  = false,          -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1,2,3},          -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   square = Plasma.Species {
      charge = 1.0, mass = mass,
      -- Velocity space grid.
      lower      = {vMin,muMin},
      upper      = {vMax,muMax},
      cells      = Nv,
      -- Initial conditions.
      init = function (t, xn)
	 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]

         return topHat(x, y, z, vpar, mu, n0, u0, vt)
      end,
      --bcx = { Plasma.Species.bcOpen,
      --        Plasma.Species.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'square'},
         frequencies = {nu},
      },
   },

   -- Neutral species with a bump in the tail.
   bump = Plasma.Species {
      charge = 1.0, mass = mass,
      -- Velocity space grid.
      lower      = {vMin,muMin},
      upper      = {vMax,muMax},
      cells      = Nv,
      -- Initial conditions.
      init = function (t, xn)
   	 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]

         return bumpMaxwell(x,y,z,vpar,mu,n0,u0,vt,ab,ub,sb,vtb)
      end,
      --bcx = { Plasma.Species.bcOpen,
      --        Plasma.Species.bcOpen },
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'bump'},
         frequencies = {nu},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve      = true, -- Evolve fields?
      externalPhi = function (t, xn) return 0.0 end,
   },
   
   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },

}
-- Run application.
plasmaApp:run()
