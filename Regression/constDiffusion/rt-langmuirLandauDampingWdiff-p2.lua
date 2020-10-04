-- Gkyl -------------------------------------------------------------------
--
-- 1x1v Landau damping simulation with diffusion in velocity space.
--
---------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Maxwellian in 1x1v, you may not need this but it made the example file more compact
local function maxwellian1D(n, vx, ux, mass, temp)
   local v2 = (vx - ux)^2
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-mass*v2/(2*temp))
end

-- Normalization parameters, shouldn't need to adjust.
epsilon0   = 1.0                         -- Permittivity of free space.
mu0        = 1.0                         -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0)   -- Speed of light.

elcMass   =  1.0   -- Electron mass.
elcCharge = -1.0   -- Electron charge.

nuElc = 0.0   --Sample electron collision frequency.

nElc = 1.0   -- Number density.
Te   = 1.0   -- Electron temperature.

-- Thermal speed, plasma frequency and Debye length.
vtElc   = math.sqrt(Te/elcMass)
wpe     = math.sqrt(elcCharge^2*nElc/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- Parameters for perturbation.
knumber      = 0.5/lambdaD
perturbation = 1e-4

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 5.0/wpe,              -- End time.
   nFrame      = 1,                    -- Number of output frames.
   lower       = {-math.pi/knumber},   -- Configuration space lower left.
   upper       = { math.pi/knumber},   -- Configuration space upper right.
   cells       = {64},                 -- Configuration space cells.
   basis       = "serendipity",        -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                    -- Polynomial order.
   cflFrac     = 0.25,
   timeStepper = "rk3",                -- One of "rk2" or "rk3".

   -- Boundary conditions for configuration space.
   periodicDirs = {1},   -- Periodic directions.
   
   -- Decomposition for configuration space.
   decompCuts = {1},     -- Cuts in each configuration direction.
   useShared  = false,   -- If to use shared memory.

   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.001,

   -- Electrons.
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- Velocity space grid.
      lower      = {-6.0*vtElc},
      upper      = { 6.0*vtElc},
      cells      = {64},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local alpha = perturbation
         local k = knumber

	 local fv = maxwellian1D(nElc, v, 0.0, elcMass, Te) 
	 return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      diff = Plasma.Diffusion {
         coefficient   = 0.0025,
         -- Optional inputs:
         diffusiveDirs = {2},
--         order         = 2,   -- Diffusion order: 2, 4, or 6.
      },
   },

   -- Fixed ions.
   ion = Plasma.Species {
      charge = -elcCharge, mass = elcMass*1823.,
      -- Velocity space grid.
      lower      = {-6.0*vtElc*math.sqrt(1./1823.)},
      upper      = { 6.0*vtElc*math.sqrt(1./1823.)},
      cells      = {64},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local fv = maxwellian1D(nElc, v, 0.0, elcMass*1823., Te) 
	 return fv
      end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k     = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
