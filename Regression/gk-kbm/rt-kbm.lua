-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"
local Logger = require "Lib.Logger"
local math = require("sci.math").generic

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2*Constants.PROTON_MASS -- (deuterium ions)
mu0 = Constants.MU0
B0 = 1.91   -- maximum field strength on axis [T]
R0 = 1.313  -- major radius of magnetic axis [m]
a  = 0.4701 -- minor radius [m]
Ti0 = 2072*eV
eps_n = 0.2
eta_e = 2.0
eta_i = 2.5
ky_rhoi = 0.5
kz_Ln = 0.1
beta_i = .02
tau = 1 -- = Ti0/Te0
-- derived parameters
r0      = 0.5*a -- minor radius of flux tube 
R       = R0 + r0 -- major radius of flux tube 
B       = B0*R0/R -- magnetic field strength in flux tube 
n0     = beta_i*B^2/(2*mu0*Ti0)
Te0     = Ti0/tau
vte  	= math.sqrt(Te0/me)
vti  	= math.sqrt(Ti0/mi)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B/mi)
omega_ce = math.abs(qe*B/me)
rho_s   = c_s/omega_ci
rho_e   = vte/omega_ce
rho_i   = vti/omega_ci
ky_min  = ky_rhoi / rho_i
dr      = 2*math.pi/ky_min
L_n     = R*eps_n
L_Te     = L_n/eta_e
L_Ti     = L_n/eta_i
kz_min  = kz_Ln / L_n
L_parallel = 2*math.pi/kz_min
-- velocity grid parameters
N_VPAR, N_MU = 8, 4
VPAR_UPPER = 4*vte
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = 16*me*vte*vte/B/2

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 3*L_n/vti, -- end time
   nFrame = 1, -- number of output frames
   lower = {r0 - .01*dr/2, -dr/2, -L_parallel/2}, -- configuration space lower left
   upper = {r0 + .01*dr/2,  dr/2, L_parallel/2}, -- configuration space upper right
   mapc2p = function(xc)
      -- field-aligned coordinates (x,y)
      local x, y, z = xc[1], xc[2], xc[3]
      -- cylindrical coordinates (R,phi)
      local R = x+R0
      local phi = y/(R0+r0)
      -- cartesian coordinates (X,Y)
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      local Z = z
      return X, Y, Z
   end,
   cells = {1, 8, 8}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2,3}, -- periodic directions

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = qe,
      mass = me,
      -- velocity space grid
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      -- initial conditions
      initBackground = Plasma.Gyrokinetic.MaxwellianProjection {
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
              isBackground = true,
             },
      init = Plasma.Gyrokinetic.MaxwellianProjection {
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
      fluctuationBCs = true, -- only apply BCs to fluctuations
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkBeta"},
   },

   -- gyrokinetic ions
   ion = Plasma.GkSpecies {
      charge = qi,
      mass = mi,
      -- velocity space grid
      lower = {VPAR_LOWER*vti/vte, MU_LOWER*mi*vti*vti/me/vte/vte},
      upper = {VPAR_UPPER*vti/vte, MU_UPPER*mi*vti*vti/me/vte/vte},
      cells = {N_VPAR, N_MU},
      -- initial conditions
      initBackground = Plasma.Gyrokinetic.MaxwellianProjection {
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
              isBackground = true,
             },
      init = Plasma.Gyrokinetic.MaxwellianProjection {
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
      fluctuationBCs = true, -- only apply BCs to fluctuations
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkBeta"},
      --diagnosticMoments = {"GkM0", "GkM1", "GkM2", perturbed = true},
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve fields?
      isElectromagnetic = true,
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return B0*R0/(R0 + x)
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()

