-- Gkyl --------------------------------------------------------------
--
-- Emission BC test.
-- This tests the external BC option with the Bronold & Fehske
-- quantum model for the electron reflection function.
--
----------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- SI units.
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q0 = 1.6021766e-19
local q_e, q_i = -q0, q0
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e17, 1.0e17
local T_e, T_i = 10*q0, q0

-- Plasma parameters.
local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local uB           = math.sqrt(T_e/m_i)
local omega_pe     = math.sqrt((n_e*q_e^2)/(epsilon_0*m_e))
local lambda_D     = math.sqrt((epsilon_0*T_e)/(n_e*q_e^2))

-- Artificially decrease the speed of light.
local c    = 10.0*vth_e
local mu_0 = 1.0/(epsilon_0*c^2)

-- Initialization function.
local function maxwellian(n, u, vth, v)
   return n/math.sqrt(2*math.pi*vth*vth)*math.exp(-(v - u)^2/(2*vth*vth))
end

sim = Plasma.App {
   logToFile = false,

   tEnd        = 10/omega_pe,       -- End time.
   nFrame      = 1,                 -- Number of output frames.
   lower       = {0.0},              -- Configuration space lower left.
   upper       = {128.0*lambda_D},   -- Configuration space upper right.
   cells       = {128},              -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                  -- Polynomial order.
   cflFrac     = 1.0,                -- CFL "fraction". Usually 1.0.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},    -- Cuts in each configuration direction.
   useShared  = false,  -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {},   -- Periodic directions.

   -- Electrons.
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- Velocity space grid.
      lower = {-6.0*vth_e},
      upper = { 6.0*vth_e},
      cells = {32},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
	 return maxwellian(n_e, 0.0, vth_e, v)
      end,
      evolve = true,   -- Evolve species?
      bcx = { Plasma.ReflectBC{},
              Plasma.BronoldFehskeBC{
                 electronAffinity = 1.0,  elemCharge   = q0,
                 effectiveMass    = 0.4,  electronMass = m_e
              } },
      diagnostics = { "M0", "M1i", "M2", "M3i", "VtSq", "Udrift", "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },

   -- Ions.
   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- Velocity space grid.
      lower = {-6.0*uB},
      upper = { 6.0*uB},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian(n_i, 0.0, vth_i, v)
      end,
      evolve = true,   -- Evolve species?
      bcx = { Plasma.ReflectBC{},
              Plasma.AbsorbBC{} },
      diagnostics = { "M0", "M1i", "M2", "M3i" },
   },
   
   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         local x =  xn[1]
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,   -- Evolve field?
      bcx = { Plasma.Field.bcSymmetry,
              Plasma.Field.bcReflect },
   },
}
-- Run application.
sim:run()
