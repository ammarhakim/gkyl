-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic()

-- This test relaxes a rectangular/square IC and a bump in tail IC.
-- Maxwellian's for comparison with each are also created.

polyOrder    = 1
nu           = 0.01                            -- Collision frequency.
n0           = 1.0                             -- Density.
u0           = 0.0                             -- Flow speed.
B0           = 1.0                             -- Background magnetic field
mass         = 1.0                             -- Define mass globally
Nx           = {2}                             -- Number of cells in configuration space.

vMin, vMax   = -2.5, 2.5                            -- Min/Max velocity in grid.
muMin, muMax = 0.0, mass*(2.5^2)/(2*B0) 
Nv           = {20,20}                         -- Number of cells in velocity space.

-- Large bump on tail of Maxwellian:
vt   = math.sqrt(1.0/12.0)                  -- Thermal speed of Maxwellian in bump.
ab   = math.sqrt(0.1)                       -- Amplitude of bump.
ub   = 4*math.sqrt(1/12.0)                    -- Location of bump.
sb   = 0.12                                 -- Softening factor to avoid divergence.
vtb  = 1.0/math.sqrt(2.0)                   -- Thermal speed of Maxwellian in bump.

-- Top hat function without drift (u=0).
local function topHat(x, vpar, mu, n, u, vth)
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
local function bumpMaxwell(x,vpar,mu,n,u,vth,bA,bU,bS,bVth)
   local Pi   = math.pi
   local vSq  = ((vpar-u)^2+2*B0*math.abs(mu)/mass)/((math.sqrt(2.0)*vth)^2)
   local vbSq = ((vpar-u)^2+2*B0*math.abs(mu)/mass)/((math.sqrt(2.0)*bVth)^2)

   return (n/math.sqrt(2.0*Pi*(vth^2))^(3/2))*math.exp(-vSq)
      +(n/math.sqrt(2.0*Pi*(bVth^2))^(3/2))*math.exp(-vbSq)*(bA^2)/((vpar-bU)^2+bS^2)
end

local maskFunc = function(t, xn)
   x, vpar, mu = xn[1], xn[2], xn[3]
   -- Mask vpar and mu space.
   local vMinMask, vMaxMask   = -2.0, 2.0     -- Min velocity in grid.
   local muMinMask, muMaxMask =  0.0, mass*(5)/(2*B0)
   ---- Mask vpar space only.
   --local vMinMask, vMaxMask   = -2.0, 2.0     -- Min velocity in grid.
   --local muMinMask, muMaxMask =  0.0, mass*(2.5^2)/(2*B0)
   ---- Mask mu space only.
   --local vMinMask, vMaxMask   = -2.5, 2.5     -- Min velocity in grid.
   --local muMinMask, muMaxMask =  0.0, mass*(5)/(2*B0)
   if (vpar < vMinMask or vpar > vMaxMask)
      or (mu < muMinMask or mu > muMaxMask) then
      return -1
   else
      return 1
   end
end

tFinal  = 10*10
nFrames = 10*1
--tFinal  = 4.146083e-02
--nFrames = 1

plasmaApp = Plasma.App {
   logToFile = false,
   
   tEnd        = tFinal,           -- End time.
   nFrame      = nFrames,             -- Number of frames to write.
   lower       = {-0.5},         -- Configuration space lower coordinate.
   upper       = { 0.5},         -- Configuration space upper coordinate.
   cells       = {Nx[1]},       -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,     -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".
   calcIntQuantEvery = 1./(10.*nFrames),
   
   -- Decomposition for configuration space.
   decompCuts = {1},            -- Cuts in each configuration direction.
   useShared  = false,          -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},          -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
--   square = Plasma.Species {
--      charge = 1.0,  mass = mass,
--      -- Velocity space grid.
--      lower = {vMin,muMin},
--      upper = {vMax,muMax},
--      cells = Nv,
--      -- Initial conditions.
--      init = function (t, xn)
--	 local x, vpar, mu = xn[1], xn[2], xn[3]
--
--         return topHat(x, vpar, mu, n0, u0, vt)
--      end,
--      -- Evolve species?
--      evolve              = true,
--      evolveCollisionless = false,
--      diagnostics         = { "M0", "M1", "M2", "intM0", "intM1", "intM2", "Upar", "VtSq" },
--      mask = maskFunc,
--      -- Collisions.
--      coll = Plasma.LBOCollisions {
--         collideWith = {'square'},
--         frequencies = {nu},
--      },
--   },

   -- Neutral species with a bump in the tail.
   bump = Plasma.Species {
      charge = 1.0,  mass = mass,
      -- Velocity space grid.
      lower = {vMin,muMin},
      upper = {vMax,muMax},
      cells = Nv,
      -- Initial conditions.
      init = function (t, xn)
   	 local x, vpar, mu = xn[1], xn[2], xn[3]

         return bumpMaxwell(x,vpar,mu,n0,u0,vt,ab,ub,sb,vtb)
      end,
      -- Evolve species?
      evolve              = true,
      evolveCollisionless = false,
      diagnostics         = { "M0", "M1", "M2", "intM0", "intM1", "intM2", "Upar", "VtSq" },
      mask = maskFunc,
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'bump'},
         frequencies = {nu},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve      = false, -- Evolve fields?
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
