-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic()

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

n0        = 1.0                             -- Density.
u0        = 0.0                             -- Flow speed.
vt        = 1.0/3.0                         -- Thermal speed..
nu_s      = 0.04                            -- Collision frequency.
nu_sb     = 0.04                            -- Collision frequency.
nu_b      = 0.02                            -- Collision frequency.
nu_bs     = 0.01                            -- Collision frequency.
B0        = 1.0                             -- Background magnetic field
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

-- Sinusoidallly perturbed Maxwellian.
local function sineMaxwell(x,vx,n,u,vth)
   local Pi   = math.pi
   local vSq  = ((vx-u)/(math.sqrt(2.0)*vth))^2
   
   return (n/math.sqrt(2.0*Pi*vth))*math.exp(-vSq)
end

plasmaApp = Plasma.App {
   logToFile = false,

   tEnd        = 100,           -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower coordinate.
   upper       = {1.0},         -- Configuration space upper coordinate.
   cells       = {8},           -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space.
   decompCuts = {1},            -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},          -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   square = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower      = {-8.0*vt},
      upper      = { 8.0*vt},
      cells      = {32},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         return topHat(x, v, n0, u0, vt)*(1.0+0.1*math.cos(2.*math.pi*x))
      end,
      evolve = true,
      diagnostics = { "M0", "M1", "M2" },
      coll = Plasma.BGKCollisions {
         collideWith = {'square', 'bump'},
         frequencies = {nu_s, nu_sb},
         exactIterFixM012 = true,
      },
   },

   -- Neutral species with a bump in the tail.
   bump = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower      = {-8.0*vt},
      upper      = { 8.0*vt},
      cells      = {32},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         return bumpMaxwell(x,v,n0,u0,vt,ab,ub,sb,vtb)*(1.0+0.1*math.cos(2.*math.pi*x))
         --return bumpMaxwell(x,v,n0,u0,vt,ab,ub,sb,vtb)
      end,
      evolve            = true,
      evolveCollisionless = false,
      diagnostics = { "M0", "M1", "M2" },
      coll = Plasma.BGKCollisions {
         collideWith = {'bump', 'square'},
         frequencies = {nu_b, nu_bs},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve      = false, -- Rvolve fields?
      externalPhi = function (t, xn) return 0.0 end,
      kperpSq     = 0.0
   },
   
   -- Magnetic geometry.
   externalField = Plasma.Geometry {
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
