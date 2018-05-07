-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- SI units
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q_e, q_i = -1.6021766e-19, 1.6021766e-19
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e18, 1.0e18
local vd_e, vd_i = 0.0, 0.0
local T_e, T_i = 10*1.6021766e-19, 10*1.6021766e-19

-- Normalized units
-- local epsilon_0, mu_0 = 1.0, 1.0
-- local q_e, q_i = -1.0, 1.0
-- local m_e, m_i = 1.0, 1836.0
-- local n_e, n_i = 1.0, 1.0
-- local vd_e, vd_i = 0.0, 0.0
-- local T_e, T_i = 1.0, 1.0

local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local omega_pe = math.sqrt((n_e * q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0 * T_e)/(n_e * q_e^2))
-- artificially decrease the speed of light
mu_0 = 1.0/(epsilon_0 * (10*vth_e)^2)

sim = Plasma.App {
   logToFile = false,

   tEnd = 100/omega_pe, -- end time
   nFrame = 10, -- number of output frames
   lower = {0.0*lambda_D}, -- configuration space lower left
   upper = {128.0*lambda_D}, -- configuration space upper right
   cells = {128}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   cflFrac = 1.0, -- CFL "fraction". Usually 1.0
   timeStepper = "rk3s4", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Plasma.VlasovSpecies {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-6.0*vth_e},
      upper = {6.0*vth_e},
      cells = {16},
      decompCuts = {1},
      -- initial conditions
      init = {"maxwellian", 
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
      evolve = true, -- evolve species?
      bcx = { Plasma.VlasovSpecies.bcAbsorb,
	      Plasma.VlasovSpecies.bcReflect },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },

   ion = Plasma.VlasovSpecies {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-6.0*vth_i},
      upper = {6.0*vth_i},
      cells = {16},
      decompCuts = {1},
      -- initial conditions
      init = {"maxwellian", 
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
      evolve = true, -- evolve species?
      bcx = { Plasma.VlasovSpecies.bcAbsorb,
	      Plasma.VlasovSpecies.bcReflect },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },
   
   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.MaxwellField.bcReflect,
	      Plasma.MaxwellField.bcSymmetry },
   },
}
-- run application
sim:run()
