-- Gkyl --------------------------------------------------------------
--
-- Basic sheath simulation using the Vlasov-Maxwell solver.
--
----------------------------------------------------------------------
local Plasma = (require "App.PlasmaOnCartGrid").VlasovMaxwell()

-- SI units.
local eV              = 1.6021766e-19
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q_e, q_i        = -1.6021766e-19, 1.6021766e-19
local m_e, m_i        = 9.109383e-31, 1.6726218e-27
local n_e, n_i        = 1.0e18, 1.0e18
local vd_e, vd_i      = 0.0, 0.0
local T_e, T_i        = 10*1.6021766e-19, 10*1.6021766e-19

-- Normalized units.
-- local epsilon_0, mu_0 = 1.0, 1.0
-- local q_e, q_i        = -1.0, 1.0
-- local m_e, m_i        = 1.0, 1836.0
-- local n_e, n_i        = 1.0, 1.0
-- local vd_e, vd_i      = 0.0, 0.0
-- local T_e, T_i        = 1.0, 1.0

local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local omega_pe     = math.sqrt((n_e * q_e^2)/(epsilon_0*m_e))
local lambda_D     = math.sqrt((epsilon_0 * T_e)/(n_e * q_e^2))
-- Artificially decrease the speed of light.
mu_0 = 1.0/(epsilon_0 * (10*vth_e)^2)

-- Collision parameters.
nuFrac = 0.1
logLambdaElc = 6.6 - 0.5*math.log(n_e/1e20) + 1.5*math.log(T_e/eV)
nuElc     = nuFrac*logLambdaElc*eV^4*n_e/(6*math.sqrt(2)*math.pi^(3/2)*epsilon_0^2*math.sqrt(m_e)*(T_e)^(3/2))  --collision freq

logLambdaIon = 6.6 - 0.5*math.log(n_i/1e20) + 1.5*math.log(T_i/eV)
nuIon     = nuFrac*logLambdaIon*eV^4*n_i/(12*math.pi^(3/2)*epsilon_0^2*math.sqrt(m_i)*(T_i)^(3/2))

nuElcIon = math.sqrt(2)*nuElc
nuIonElc = nuElcIon/(m_i/m_e)

sim = Plasma.App {
   logToFile = false,

   tEnd        = 100/omega_pe,     -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0.0*lambda_D},   -- Configuration space lower left.
   upper       = {128.0*lambda_D}, -- Configuration space upper right.
   cells       = {64},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3s4",          -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   parallelizeSpecies = true,

   -- Boundary conditions for configuration space.
   periodicDirs = {}, -- Periodic directions.

   -- Electrons.
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- Velocity space grid.
      lower = {-6.0*vth_e},
      upper = { 6.0*vth_e},
      cells = {16},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            return n_e
         end,
         temperature = function (t, xn)
            return T_e
         end,
         driftSpeed = function (t, xn)
            return vd_e
         end
      },
--      coll = Plasma.LBOCollisions {
--        collideWith = {'elc', 'ion'},
--        normNu      = {nuElc*((2*(vth_e^2))^(3/2))/n_e, nuElcIon*((vth_e^2+vth_i^2)^(3/2))/n_e},
--      },
      evolve = true, -- Evolve species?
      bcx = { Plasma.AbsorbBC{diagnostics={"intM0", "intM1i", "intM2"}},
              Plasma.ReflectBC{diagnostics={"intM0", "intM1i", "intM2"}} },
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
   },

   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- Velocity space grid.
      lower = {-6.0*vth_i},
      upper = { 6.0*vth_i},
      cells = {16},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            return n_i
         end,
         temperature = function (t, xn)
            return T_i
         end,
         driftSpeed = function (t, xn)
            return vd_i
         end
      },
--      coll = Plasma.LBOCollisions {
--        collideWith = {'ion', 'elc'},
--        normNu      = {nuIon*((2*(vth_i^2))^(3/2))/n_i, nuIonElc*((vth_i^2+vth_e^2)^(3/2))/n_i},
--      },
      evolve = true, -- Evolve species?
      bcx = { Plasma.AbsorbBC{diagnostics={"intM0", "intM1i", "intM2"}},
              Plasma.ReflectBC{diagnostics={"intM0", "intM1i", "intM2"}} },
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
   },
   
   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
      bcx = { Plasma.Field.bcReflect,
              Plasma.Field.bcSymmetry },
   },
}
-- Run application.
sim:run()
