-- Gkyl -------------------------------------------------------------------
--
-- GK ion sound wave test case
--
---------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").Gyrokinetic()

TiTe    = 1.0
knumber = 0.5
XL, XU  = -math.pi/knumber, math.pi/knumber
Ti0     = 1.0
Te0     = Ti0/TiTe
n0      = 1.0
ni0     = n0
ne0     = n0
B0      = 1.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 15,                 -- End time.
   nFrame      = 1,                  -- Number of output frames.
   lower       = {-math.pi/knumber}, -- Configuration space lower left.
   upper       = {math.pi/knumber},  -- Configuration space upper right.
   cells       = {16},               -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".
   cflFrac     = 1.0,

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = 1.0,  mass = 1.0,
      -- Velocity space grid
      lower = {-6.0},
      upper = { 6.0},
      cells = {32},
      -- Initial conditions.
      -- Specify background so that we can plot perturbed distribution and moments.
      background = Plasma.MaxwellianProjection {
         density     = function (t, xn) return ni0 end,
         temperature = function (t, xn) return Ti0 end,
      },
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, v = xn[1], xn[2]
            local k = knumber
            local alpha = 0.01
            local perturb = alpha*math.cos(k*x)
            return ni0*(1+perturb)
         end,
         temperature = function (t, xn) return Ti0 end,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "M2", "intM0", "intM2"},
   },

   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = -1.0,
      mass = 1.0,
      temp = Te0,
      -- Initial conditions.. use ion background so that background is exactly neutral.
      init = function (t, xn)
         return ne0
      end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
      kperpSq = 0.0,
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         return B0
      end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
