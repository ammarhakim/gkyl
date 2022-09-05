-- Gkyl -------------------------------------------------------------------
--
-- GK ion acoustic wave damping.
-- Time normalized to the ion gyrofrequency, omega_ci 
-- Distances normalized to the ion sound gyroradius, rho_s=c_s/omega_ci.
-- Temperatures normalized to the electron temperature, T_e.
-- Velocities normalized to c_s=sqrt(T_e/m_i).
--
---------------------------------------------------------------------------

local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic()

tau     = 1.0            -- Ion temperature to electron temperature ratio.
Te      = 1.0            -- Electron temperature.
mIon    = 1.0            -- Ion mass
mElc    = 1.0/1836.16    -- Electron mass.
n0      = 1.0            -- Background density.
B0      = 1.0            -- Background magnetic field magnitude.
nuIon   = 2.0            -- Ion-ion collision frequency.
knumber = 0.5            -- Wave number.

-- Lower and upper bound of configuration space.
xLower, xUpper = -math.pi/knumber, math.pi/knumber

Ti      = Te*tau     -- Ion temperature.
nIon    = n0         -- Background ion density.
nElc    = n0         -- Background electron density.

vtIon   = math.sqrt(Ti)                -- Ion thermal speed. 
vtElc   = math.sqrt((mIon/mElc)*Te)    -- Electron thermal speed. 

tFinal  = 15        -- Final simulation time.

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = tFinal,           -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {xLower},         -- Configuration space lower left.
   upper       = {xUpper},         -- Configuration space upper right.
   cells       = {8},              -- Configuration space cells.
   mapc2p = function (xc)
      return xc[1]
   end,
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 1.0,

   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.001,

   -- Decomposition for configuration space.
   decompCuts = {1},    -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},   -- periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = 1.0,  mass   = mIon,
      -- Velocity space grid.
      lower      = {-6.0*vtIon, 0.0},
      upper      = { 6.0*vtIon, mIon*((5.0*vtIon)^2)/(2.0*B0)},
      cells      = {64, 8},
      decompCuts = {1, 1},
      -- Initial conditions.
      -- Specify background so that we can plot perturbed distribution and moments.
      background = Plasma.MaxwellianProjection{
         density     = function (t, xn) return nIon end,
         temperature = function (t, xn) return Ti end,
      },
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x, v, mu = xn[1], xn[2], xn[3]
            local k        = knumber
            local alpha    = 0.01
            local perturb  = alpha*math.cos(k*x)
            return nIon*(1.0+perturb)
         end,
         temperature = function (t, xn) return Ti end,
      },
--      coll = Plasma.LBOCollisions{
--         collideWith = { 'ion' },
--         frequencies = { nuIon },
--      }, 
      polarizationDensityFactor = nIon,
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "M2"},
   },

   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = -1.0,
      mass   = mElc,
      temp   = Te,
      -- Initial conditions.. use ion background so that background is exactly neutral.
      init = function (t, xn) return nElc end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
      kperpSq = 0.0,
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
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
