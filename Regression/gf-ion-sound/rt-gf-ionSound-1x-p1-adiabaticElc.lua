-- Gkyl -------------------------------------------------------------------
--
-- Gyrofluid ion acoustic wave damping with adiabatic electrons.
-- Time normalized to the ion gyrofrequency, omega_ci 
-- Distances normalized to the ion sound gyroradius, rho_s=c_s/omega_ci.
-- Temperatures normalized to the electron temperature, T_e.
-- Velocities normalized to c_s=sqrt(T_e/m_i).
--
---------------------------------------------------------------------------

local Plasma = require("App.PlasmaOnCartGrid").Gyrofluid()

tau     = 1.0            -- Ion temperature to electron temperature ratio.
Te      = 1.0            -- Electron temperature.
mIon    = 1.0            -- Ion mass
mElc    = 1.0/1836.16    -- Electron mass.
n0      = 1.0            -- Background density.
B0      = 1.0            -- Background magnetic field magnitude.
kpar0   = 0.5            -- Wave number.
kperp   = 0.3            -- k_perp*rho_s.
nuIon   = 0.0            -- Ion-ion collision frequency.
nuElc   = 0.0            -- Electron-electron collision frequency.

-- Lower and upper bound of configuration space.
xLower, xUpper = -math.pi/kpar0, math.pi/kpar0

Ti      = Te*tau     -- Ion temperature.
nIon    = n0         -- Background ion density.
nElc    = n0         -- Background electron density.

vtIon   = math.sqrt(Ti)                -- Ion thermal speed. 
vtElc   = math.sqrt((mIon/mElc)*Te)    -- Electron thermal speed. 

-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

kappaParIon  = nIon*(vtIon^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vtIon*kpar0+nuIon)
kappaPerpIon = nIon*(vtIon^2)/(math.sqrt(2)*D_perp*vtIon*kpar0+nuIon)
kappaParElc  = nElc*(vtElc^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vtElc*kpar0+nuElc)
kappaPerpElc = nElc*(vtElc^2)/(math.sqrt(2)*D_perp*vtElc*kpar0+nuElc)

plasmaApp = Plasma.App {
   tEnd        = 5.0,              -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {xLower},         -- Configuration space lower left.
   upper       = {xUpper},         -- Configuration space upper right.
   cells       = {16},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 0.90,

   -- Frequency with which to compute integrated diagnostics.
--   calcIntQuantEvery = 1/600.,

   -- Decomposition for configuration space.
   decompCuts = {1},    -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},   -- Periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = 1.0,  mass = mIon,
      -- Initial conditions.
      -- Specify background so that we can plot perturbed distribution and moments.
      background = Plasma.GyrofluidProjection{
         density = function (t, xn) return nIon end,
         parallelTemperature = function (t, xn) return Ti end,
         perpendicularTemperature = function (t, xn) return Ti end,
      },
      init = Plasma.GyrofluidProjection{
         density = function (t, xn)
            local x     = xn[1]
            local k     = kpar0
            local alpha = 0.01
            return nIon*(1.0+alpha*math.cos(k*x))
         end,
         parallelTemperature = function (t, xn) return Ti end,
         perpendicularTemperature = function (t, xn) return Ti end,
      },
      closure = Plasma.HeatFlux{
         kappaPar = kappaParIon,  kappaPerp = kappaPerpIon,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"M0","intM0","intM1","intM2","Upar","Tpar","Tperp"},
   },

   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = -1.0,  mass = mElc,
      temp   = Te,
      -- Initial conditions.. use ion background so that background is exactly neutral.
      init = function (t, xn) return nElc end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
      kperpSq = kperp^2,
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn) return B0 end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
