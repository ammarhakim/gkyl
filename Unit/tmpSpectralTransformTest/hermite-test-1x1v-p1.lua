-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- This test initializes a Maxwellian weighted by the 10th Hermite polynomial
-- For testing Hermite transform post-processing and eigenfunctions of the collision operator

polyOrder = 1       -- Polynomial order.
n0        = 1.0     -- Plasma density.
charge    = 0.0     -- Species charge.
u0        = 0.0     -- Flow speed.
vt        = 1.0     -- Thermal speed.
mass      = 1.0     -- Species mass.
nu        = 0.01    -- Collision frequency.

local function factorial(n)
   if (n==0) then
      return 1
   else
      return n*factorial(n-1)
   end
end

-- Maxwellian in 1x1v weighted by 5th Hermite polynomial
local function maxwellian1DHermite5(n, vx, ux, temp, mass)
   local vnorm = (vx - ux)/math.sqrt(2.0*temp/mass)
   local hermite5 = (32*vnorm^5 - 160*vnorm^3 + 120*vnorm)/math.sqrt(math.pow(2, 5)*factorial(5))
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-vnorm^2)*(1+hermite5)
end

-- Maxwellian in 1x1v weighted by 10th Hermite polynomial
local function maxwellian1DHermite10(n, vx, ux, temp, mass)
   local vnorm = (vx - ux)/math.sqrt(2.0*temp/mass)
   local hermite10 = (1024*vnorm^10 - 23040*vnorm^8 + 161280*vnorm^6 - 403200*vnorm^4 + 302400*vnorm^2 - 30240)/math.sqrt(math.pow(2, 10)*factorial(10))
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-vnorm^2)*(1+hermite10)
end

-- Maxwellian in 1x1v weighted by 5th and 10th Hermite polynomial
local function maxwellian1DMultHermite(n, vx, ux, temp, mass)
   local vnorm = (vx - ux)/math.sqrt(2.0*temp/mass)
   local hermite5 = (32*vnorm^5 - 160*vnorm^3 + 120*vnorm)/math.sqrt(math.pow(2, 5)*factorial(5))
   local hermite10 = (1024*vnorm^10 - 23040*vnorm^8 + 161280*vnorm^6 - 403200*vnorm^4 + 302400*vnorm^2 - 30240)/math.sqrt(math.pow(2, 10)*factorial(10))
   local hermite20 = (1048576*vnorm^20 - 99614720*vnorm^18 + 3810263040*vnorm^16 - 76205260800*vnorm^14 + 866834841600*vnorm^12 - 5721109954560*vnorm^10 + 21454162329600*vnorm^8 - 42908324659200*vnorm^6 + 40226554368000*vnorm^4 - 13408851456000*vnorm^2 + 670442572800)/math.sqrt(math.pow(2, 20)*factorial(20))
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-vnorm^2)*(1+hermite20+hermite10+hermite5)
end

local function hermiteBasis(mmax, x)
  local hmx = { [0] = math.exp(-(x^2))/math.sqrt(math.pi),
                      math.exp(-(x^2))*math.sqrt(2.0/math.pi)*x
  }

  local function _hermiteBasis(m,x)
     if m == 0 then
        return hmx[0]
     elseif m == 1 then
        return hmx[1]
     end
     return (math.sqrt(2.0/m))*x*hmx[m-1]-math.sqrt((m-1)/m)*hmx[m-2]
  end

  for m = 2, mmax do
     hmx[m] = _hermiteBasis(m,x)
  end
  return hmx[mmax]
end
local function maxwellian1DMultHermiteP(n, vx, ux, temp, mass)
   local hermite0  = hermiteBasis(0,vx)
   local hermite5  = hermiteBasis(5,vx)
   local hermite10 = hermiteBasis(10,vx)
   local hermite20 = hermiteBasis(20,vx)
--   return n/math.sqrt(2*math.pi*temp/mass)*(hermite0+hermite5+hermite10+hermite20)
   return (hermite0+hermite5+hermite10+hermite20)
end

plasmaApp = Plasma.App {
   logToFile = false,

   tEnd        = 100,           -- End time.
   nFrame      = 100,           -- Number of frames to write.
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
   neut = Plasma.Species {
      charge = charge, mass = mass,
      -- Velocity space grid.
      lower      = {-12.0*vt},
      upper      = {12.0*vt},
      cells      = {96},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]

         --return maxwellian1DMultHermite(n0, v, u0, vt, mass)
         return maxwellian1DMultHermiteP(n0, v, u0, vt, mass)
      end,
      -- Evolve species?
      evolve = true,
      -- Diagnostic moments.
      diagnosticMoments     = { "M0", "M1i", "M2", "vtSq" },
      phaseSpaceDiagnostics = { "Hermite2" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {'neut'},
         frequencies = {nu},
      },
   },

}
-- run application
plasmaApp:run()
