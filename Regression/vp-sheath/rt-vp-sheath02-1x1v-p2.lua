-- Gkyl --------------------------------------------------------------
--
-- Basic sheath simulation using Vlasov-Poisson.
--
----------------------------------------------------------------------
local Plasma = (require "App.PlasmaOnCartGrid").VlasovMaxwell()

-- SI units.
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
            return {vd_e}
         end
      },
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
            return {vd_i}
         end
      },
      evolve = true, -- Evolve species?
      bcx = { Plasma.AbsorbBC{diagnostics={"intM0", "intM1i", "intM2"}},
	      Plasma.ReflectBC{diagnostics={"intM0", "intM1i", "intM2"}} },
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
   },
   
   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      evolve   = true, -- Evolve field?
      hasMagneticField = false,
      bcLowerPhi = {{T="D", V=0.0}}, 
      bcUpperPhi = {{T="N", V=0.0}}, 
   },
}
-- Run application.
sim:run()
