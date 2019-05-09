-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- Might help to think of species 1 as the ions.
local n1, u1, p1, m1 = 1.0, 0.10, 1.0, 1.0
local n2, u2, p2, m2 = 1.0, 2.50, 0.50, 0.05

local vt1 = math.sqrt(p1/(n1*m1))
local vt2 = math.sqrt(p2/(n2*m2))

local K1 = 0.01/n1
local K2 = 0.01/n2

-- Collisionalities
nu11 = vt1/K1 
nu22 = vt2/K2 
nu21 = 2*vt2/K2
nu12 = (m2/m1)*nu21

print(" ")
print("tEnd = ", 0.5)
print(" ")
print('1-1 collision period: ', 1.0/nu11)
print('2-2 collision period: ', 1.0/nu22)
print('1-2 collision period: ', 1.0/nu12)
print('2-1 collision period: ', 1.0/nu21)
print(' ')

vlasovApp = Plasma.App {
   logToFile = false,

   tEnd        = 0.005,            -- End time.
   nFrame      = 1,                -- Number of frames to write.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {1.0},            -- Configuration space upper right.
   cells       = {16},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".
--   cflFrac     = 1.0,
--   maximumDt   = 1.30e-05,

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- First neutral species.
   neut1 = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0, mass = m1,
      -- Velocity space grid.
      lower = {-6.0*vt1},
      upper = { 6.0*vt1},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return n1
         end,
         driftSpeed = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return {u1}
         end,
         temperature = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return m1*(vt1^2)
         end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true,
      --evolveCollisions = false,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2", "u", "vtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      -- Collisions.
--      coll = Plasma.LBOCollisions {
--         collideWith  = { "neut1", "neut2", },
--         frequencies  = { nu11, nu12, },
----         collideWith  = { "neut1", },
----         frequencies  = { nu11, },
----         collideWith  = { "neut2" },
----         -- Optional arguments:
----         betaGreene  = 0.0,    -- Free parameter, must be >-1.
--      },
   },

   -- Second neutral species.
   neut2 = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0, mass = m2,
      -- Velocity space grid.
      lower = {-6.0*vt2},
      upper = { 6.0*vt2},
      cells = {32},
      decompCuts = {1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return n2
         end,
         driftSpeed = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return {u2}
         end,
         temperature = function (t, zn)
            local x, vpar = zn[1], zn[2]
            return m2*(vt2^2)
         end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true,
      --evolveCollisions = false,
      -- Diagnostic moments.
      diagnosticMoments = { "M0", "M1i", "M2", "u", "vtSq" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = { "neut1", "neut2", },
         frequencies = { nu21, nu22, },
--         collideWith = { "neut2", },
--         frequencies = { nu22, },
--         collideWith = { "neut1" },
--         -- Optional arguments:
--         betaGreene  = 0.0,    -- Free parameter, must be >-1.
      },
   },
}
-- Run application.
vlasovApp:run()
