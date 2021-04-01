-- Gkyl --------------------------------------------------------------
--
-- BGK Relexation test
--
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").Gyrokinetic()

n0        = 1.0                             -- Density.
B0        = 1.0                             -- Background magnetic field

sim = Plasma.App {
   logToFile = false,

   tEnd        = 10.0,             -- End time.
   nFrame      = 1,                -- Number of frames to write.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {2*math.pi},      -- Configuration space upper right.
   cells       = {4},              -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".
--   cflFrac     = 0.1,

   -- Decomposition for configuration space.
   decompCuts = {1},     -- Cuts in each configuration direction.
   useShared  = false,   -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},    -- Periodic directions.

   -- Neutral species (the charge=1 doesn't matter here).
   spec = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = { 6.0},
      cells = {16},
      -- Initial conditions.
      init = function (t, xn)
         local x, v = xn[1], xn[2]
         if math.abs(v) > 1.5 then return 0 else return 1.0/3 end
      end,
      -- Evolve species?
      evolve              = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnosticMoments = { "GkM0", "GkM1", "GkM2" },
      diagnosticIntegratedMoments = { "intM0", "intM1", "intM2" },
      -- Collisions.
      bgk = Plasma.BGKCollisions {
         collideWith = {"spec"},
         frequencies = {1.0},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve      = false, -- Rvolve fields?
--      externalPhi = function (t, xn) return 0.0 end,
      kperp2      = 0.0
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
sim:run()
