-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"

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
ky_rhos = 0.01
kz_Ln = 1.0
beta_hat = 10.0 -- beta_e/2*m_i/m_e
beta_e = beta_hat*2*me/mi
tau = 1 -- Ti0/Te0
-- derived parameters
r0      = 0.5*a -- minor radius of flux tube 
R       = R0 + r0 -- major radius of flux tube 
B       = B0*R0/R -- magnetic field strength in flux tube 
Te0     = Ti0/tau
n0     = beta_e*B^2/(2*mu0*Te0)
ne0 = n0
ni0 = n0
vte  	= math.sqrt(Te0/me)
vti  	= math.sqrt(Ti0/mi)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B/mi)
omega_ce = math.abs(qe*B/me)
rho_s   = c_s/omega_ci
rho_e   = vte/omega_ce
rho_i   = vti/omega_ci
ky_min  = ky_rhos / rho_s
dr      = 2*math.pi/ky_min
L_n     = R*eps_n
L_Te     = L_n/eta_e
L_Ti     = L_n/eta_i
kz_min  = kz_Ln / L_n
L_parallel = 2*math.pi/kz_min
-- velocity grid parameters
VPAR_UPPER = 6*vte
VPAR_LOWER = -VPAR_UPPER
omega = math.sqrt(kz_Ln^2/beta_hat/(1+ky_rhos^2/beta_hat))*vte/L_n
omega = math.sqrt(kz_min^2*2/beta_e*c_s^2/(1+beta_hat*ky_rhos^2))

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 2*math.pi/omega, --5/knumber, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/kz_min}, -- configuration space lower left
   upper = {math.pi/kz_min}, -- configuration space upper right
   cells = {16}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = qe,
      mass = me,
      -- velocity space grid
      lower = {-6*vte},
      upper = {6*vte},
      cells = {32},
      decompCuts = {1},
      -- initial conditions
      initBackground = Plasma.Gyrokinetic.MaxwellianProjection {
         density = function (t, xn)
            local x = xn[1]
            return ne0
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0
         end,
         exactScaleM0 = true,
         isBackground = true,
      },
      init = Plasma.Gyrokinetic.MaxwellianProjection {
         density = function (t, xn)
	    local x = xn[1]
            return ne0*(1 + 1e-6*math.cos(kz_min*x))
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0
         end,
         exactScaleM0 = true,
      },
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkM1", perturbed=true},
   },

   ion = Plasma.GkSpecies {
      charge = qi,
      mass = mi,
      -- velocity space grid
      lower = {-6*vti},
      upper = {6*vti},
      cells = {32},
      decompCuts = {1},
      init = Plasma.Gyrokinetic.MaxwellianProjection {
         density = function (t, xn)
	    local x = xn[1]
            return ni0
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Ti0
         end,
         exactScaleM0 = true,
      },
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkM1"},
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve field?
      isElectromagnetic = true,
      kperp2 = ky_min^2,
      mu0 = mu0,
      discontinuousPhi = false,
      discontinuousApar = true,
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         return B
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
