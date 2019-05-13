-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- Normalization parameters.
epsilon0   = 1.0                       -- Permittivity of free space.
mu0        = 1.0                       -- Permiability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- Speed of light.

elcMass   =  1.0     -- Electron mass.
elcCharge = -1.0     -- Electron charge.

ionMass   = 1836.153 -- Ion mass.
ionCharge = 1.0      -- Ion charge.

nuElc    = 0.001                    -- Electron collision frequency.
nuIon    = nuElc/math.sqrt(ionMass) -- Ion collision frequency.
nuElcIon = nuElc/1.96               -- Electron-ion collision frequency.
nuIonElc = elcMass*nuElcIon/ionMass -- Ion-electron collision frequency.

Te_Ti = 1.0    -- Ratio of electron to ion temperature.

n  = 1.0
Te = 1.0
Ti = 1.0/Te_Ti

-- Derived parameters.
vtElc = math.sqrt(Te/elcMass)
vtIon = math.sqrt(Ti/ionMass)
-- Plasma frequency and Debye length.
wpe     = math.sqrt((elcCharge^2)*n/(epsilon0*elcMass))
wpi     = math.sqrt((ionCharge^2)*n/(epsilon0*ionMass))
lambdaD = vtElc/wpe

-- Domain size and simulation time.
LX = 64.0*lambdaD

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd         = 1.0/nuElcIon,  -- End time.
   nFrame       = 1,             -- Number of output frames.
   lower        = {-LX/2.0},     -- Configuration space lower left.
   upper        = { LX/2.0},     -- Configuration space upper right.
   cells        = {8},           -- Configuration space cells.
   basis        = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder    = 1,             -- Polynomial order.
   timeStepper  = "rk3",         -- One of "rk2" or "rk3".
   -- Boundary conditions for configuration space.
   periodicDirs = {1},           -- Periodic directions.
   -- Decomposition for configuration space.
   decompCuts   = {1},           -- Cuts in each configuration direction.
   useShared    = false,         -- If to use shared memory.

   -- Integrated moment flag, compute quantities 1000 times in simulation.
--   calcIntQuantEvery = 0.0025,

   -- Electrons.
   elc = Plasma.Species {
      charge = 0.0, mass = elcMass,
      -- Velocity space grid.
      lower = {-5.0*vtElc, -5.0*vtElc},
      upper = { 5.0*vtElc,  5.0*vtElc},
      cells = {12, 12},
      decompCuts = {1, 1}, -- Do not change, no parallelization in velocity space currently.
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            local xm = -0.125*LX
            if x>xm then
               return n*(1+math.exp(-0.075*(x-xm)^2))
            else
               return n*(1+math.exp(-0.75*(x-xm)^2))
            end
         end,
         driftSpeed = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            return {0.0, 0.0}
         end,
         temperature = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            return Te
         end,
         exactScaleM0    = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?
      evolveCollisionless = false,
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i", "u", "vtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      coll = Plasma.LBOCollisions {
         collideWith  = { "elc",  "ion" },
         frequencies  = { nuElc,  nuElcIon },
--         collideWith  = { "elc" },
--         frequencies  = { nuElc },
--         collideWith  = { "ion" },
--         frequencies  = { nuElcIon },
         -- Optional parameters:
--         betaGreene   = 1.0,    -- Free parameter, must be >-1.
      },
   },
   -- Protons.
   ion = Plasma.Species {
      charge = 0.0, mass = ionMass,
      -- Velocity space grid.
      lower = {-5.0*vtIon, -5.0*vtIon},
      upper = { 5.0*vtIon,  5.0*vtIon},
      cells = {12, 12},
      decompCuts = {1, 1}, -- Do not change, no parallelization in velocity space currently.
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            local xm = -0.125*LX
            if x>xm then
               return n*(1+math.exp(-0.075*(x-xm)^2))
            else
               return n*(1+math.exp(-0.75*(x-xm)^2))
            end
         end,
         driftSpeed = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            return {vtIon, 0.0}
         end,
         temperature = function (t, zn)
            local x, vx, vy = zn[1], zn[2], zn[3]
            return Ti
         end,
         exactScaleM0    = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- Evolve species?
      evolveCollisionless = false,
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i", "u", "vtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      coll = Plasma.LBOCollisions {
         collideWith  = { "ion", "elc" },
         frequencies  = { nuIon, nuIonElc },
--         collideWith  = { "ion" },
--         frequencies  = { nuIon },
--         collideWith  = { "elc" },
--         frequencies  = { nuIonElc },
         -- Optional parameters:
--         betaGreene   = 1.0,    -- Free parameter, must be >-1.
      },
   },

}
-- Run application.
plasmaApp:run()
