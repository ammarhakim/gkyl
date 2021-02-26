-- Gkyl ------------------------------------------------------------------------
local Plasma    = require ("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Logger    = require "Lib.Logger"
local math      = require("sci.math").generic

-- Physical parameters.
eV      = Constants.ELEMENTARY_CHARGE
qe      = -eV
qi      = eV
me      = Constants.ELECTRON_MASS
mi      = 2*Constants.PROTON_MASS   -- Deuterium ions.
mu0     = Constants.MU0
B0      = 1.91   -- Maximum field strength on axis [T].
R0      = 1.313  -- Major radius of magnetic axis [m].
a       = 0.4701 -- Minor radius [m].
Ti0     = 2072*eV
eps_n   = 0.2
eta_e   = 2.0
eta_i   = 2.5
ky_rhoi = 0.5
kz_Ln   = 0.1
beta_i  = .02
tau     = 1 -- = Ti0/Te0
-- Derived parameters.
r0         = 0.5*a     -- Minor radius of flux tube.
R          = R0 + r0   -- Major radius of flux tube.
B          = B0*R0/R   -- Magnetic field strength in flux tube.
n0         = beta_i*B^2/(2*mu0*Ti0)
Te0        = Ti0/tau
vte        = math.sqrt(Te0/me)
vti        = math.sqrt(Ti0/mi)
c_s        = math.sqrt(Te0/mi)
omega_ci   = math.abs(qi*B/mi)
omega_ce   = math.abs(qe*B/me)
rho_s      = c_s/omega_ci
rho_e      = vte/omega_ce
rho_i      = vti/omega_ci
ky_min     = ky_rhoi / rho_i
dr         = 2*math.pi/ky_min
L_n        = R*eps_n
L_Te       = L_n/eta_e
L_Ti       = L_n/eta_i
kz_min     = kz_Ln / L_n
L_parallel = 2*math.pi/kz_min
-- Velocity grid parameters.
N_VPAR, N_MU = 8, 4
VPAR_UPPER   = 4*vte
VPAR_LOWER   = -VPAR_UPPER
MU_LOWER     = 0
MU_UPPER     = 16*me*vte*vte/B/2

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = 3*L_n/vti/2,   -- End time.
   nFrame = 1,           -- Number of output frames.
   lower  = {r0 - .01*dr/2, -dr/2, -L_parallel/2}, -- Configuration space lower left.
   upper  = {r0 + .01*dr/2,  dr/2, L_parallel/2},  -- Configuration space upper right.
   mapc2p = function(xc)
      -- Field-aligned coordinates (x,y).
      local x, y, z = xc[1], xc[2], xc[3]
      -- Cylindrical coordinates (R,phi).
      local R = x+R0
      local phi = z/(R0+r0)
      -- Cartesian coordinates (X,Y).
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      local Z = y
      return X, Y, Z
   end,
   cells       = {1, 8, 8},       -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2" or "rk3".
   cflFrac     = 1.0,

   -- Decomposition for configuration space.
   useShared = false,   -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1,2,3},   -- Periodic directions.

   -- Gyrokinetic electrons.
   electron = Plasma.Species {
      charge = qe,
      mass = me,
      -- Velocity space grid.
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      -- Initial conditions.
      background = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x = xn[1]
            return n0*(1-(x-r0)/L_n)
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_Te)
         end,
         exactScaleM012 = true,
      },
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            local perturb = 1e-5*rho_s/L_n*math.cos(ky_min*y+kz_min*z)
            return n0*(1-(x-r0)/L_n) + n0*perturb
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_Te)
         end,
         exactScaleM012 = true,
      },
      fluctuationBCs = true,   -- Only apply BCs to fluctuations.
      evolve         = true,   -- Evolve species?
      diagnosticMoments = {"GkM0", "GkBeta"},
      diagnosticIntegratedMoments = {"intM0", "intM2"},
   },

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = qi,
      mass = mi,
      -- Velocity space grid.
      lower = {VPAR_LOWER*vti/vte, MU_LOWER*mi*vti*vti/me/vte/vte},
      upper = {VPAR_UPPER*vti/vte, MU_UPPER*mi*vti*vti/me/vte/vte},
      cells = {N_VPAR, N_MU},
      -- Initial conditions.
      background = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x = xn[1]
            return n0*(1-(x-r0)/L_n)
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Ti0*(1-(x-r0)/L_Ti)
         end,
         exactScaleM012 = true,
      },
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            local perturb = 1e-5*rho_s/L_n*math.cos(ky_min*y+kz_min*z)
            return n0*(1-(x-r0)/L_n) --+ n0*perturb
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Ti0*(1-(x-r0)/L_Ti)
         end,
         exactScaleM012 = true,
      },
      fluctuationBCs = true,   -- Only apply BCs to fluctuations.
      evolve         = true,   -- Evolve species?
      diagnosticMoments = {"GkM0", "GkBeta"},
      diagnosticIntegratedMoments = {"intM0", "intM2"},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve            = true, -- Evolve fields?
      isElectromagnetic = true,
   },

   -- Magnetic geometry. 
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return B0*R0/(R0 + x)
      end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()

