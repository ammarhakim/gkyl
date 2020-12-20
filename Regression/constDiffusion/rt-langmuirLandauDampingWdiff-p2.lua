-- Gkyl -------------------------------------------------------------------
--
-- 1x1v Landau damping simulation with diffusion in velocity space.
--
---------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Maxwellian in 1x1v, you may not need this but it made the example file more compact
local function maxwellian1D(n, vx, ux, vt)
   local v2 = (vx - ux)^2
   return n/math.sqrt(2*math.pi*vt^2)*math.exp(-v2/(2*vt^2))
end

-- Normalization parameters, shouldn't need to adjust.
epsilon0   = 1.0                         -- Permittivity of free space.
mu0        = 1.0                         -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0)   -- Speed of light.

qe = -1.0
qi = 1.0
me = 1.0
mi = 1836.153
vte = 1.0 -- electron thermal velocity
vti = vte/math.sqrt(mi) -- ion thermal velocity 

nElc = 1.0   -- Number density.

-- plasma frequency and Debye length.
wpe = math.sqrt(qe^2*nElc/(epsilon0*me))
lambdaD = vte/wpe

-- Parameters for perturbation.
knumber      = 0.5/lambdaD
perturbation = 1.0e-2 -- distribution function perturbation
Nx = 32
Nvx = 32
vMaxElc = 6.0*vte
vMaxIon = 6.0*vti
-- Compute hyperdiffusion coefficients from dv (so coefficient is a frequency)
nuElc = 1.0e-1*(vMaxElc/Nvx)^6

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 100.0/wpe,              -- End time.
   nFrame      = 1,                    -- Number of output frames.
   lower       = {-math.pi/knumber},   -- Configuration space lower left.
   upper       = { math.pi/knumber},   -- Configuration space upper right.
   cells       = {Nx},                 -- Configuration space cells.
   basis       = "serendipity",        -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                    -- Polynomial order.
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
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {-vMaxElc},
      upper = {vMaxElc},
      cells = {Nvx},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         local alpha = perturbation
         local k = knumber

         local fv = maxwellian1D(nElc, v, 0.0, vte) 
         return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      diff = Plasma.Diffusion {
         coefficient   = nuElc,
         -- Optional inputs:
         diffusiveDirs = {2},
         order = 6,   -- Diffusion order: 2, 4, or 6.
      },
   },

   -- Fixed ions.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {-vMaxIon},
      upper = {vMaxIon},
      cells = {Nvx},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         local fv = maxwellian1D(nElc, v, 0.0, vti) 
         return fv
      end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
