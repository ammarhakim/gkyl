-- 4D GK ETG linear instability calculation
-- using 'pgkyl -f etg4d_elecEnergy_ growth' should approximately give growth rate printed at end of run
--
-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"
local math = require("sci.math").generic

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0 = 2072*eV 
Ti0 = 2072*eV 
B0 = 1.91   -- magnetic field strength on axis [T]
R0 = 1.313  -- major radius of magnetic axis [m]
a  = 0.4701 -- minor radius [m]
n0 = 4.992e19 -- [1/m^3]
-- derived parameters
r0       = 0.5*a -- minor radius of flux tube 
R        = R0 + r0 -- major radius of flux tube 
B        = B0*R0/R -- magnetic field strength in flux tube 
vte  	 = math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B/mi)
omega_ce = math.abs(qe*B/me)
rho_s    = c_s/omega_ci
rho_e    = vte/omega_ce
dr       = 32*rho_e
L_T      = R/10 
ky_min   = 2*math.pi/dr
omegade  = ky_min*rho_e*vte/R
-- velocity grid parameters
N_VPAR, N_MU = 16, 8
VPAR_UPPER = 4*vte
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = 16*me*vte*vte/B/2

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = .5e-6, -- end time
   nFrame = 1, -- number of output frames
   lower = {r0 - 0.001*dr/2, -dr/2}, -- configuration space lower left
   upper = {r0 + 0.001*dr/2,  dr/2}, -- configuration space upper right
   cells = {1, 8}, -- configuration space cells
   mapc2p = function(xc)
      -- field-aligned coordinates (x,y)
      local x, y = xc[1], xc[2]
      -- cylindrical coordinates (R,phi)
      local R = x+R0
      local phi = y/(R0+x)
      -- cartesian coordinates (X,Y)
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      return X, Y
   end,
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 0.9,

   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2}, -- periodic directions
   deltaF = true, -- only apply BCs to fluctuations, and use perturbed moments in field solve

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = qe,
      mass = me,
      -- velocity space grid
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      decompCuts = {1, 1},
      -- initial conditions
      initBackground = Plasma.Gyrokinetic.MaxwellianProjection {
         density = function (t, xn)
            local x = xn[1]
            return n0
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_T)
         end,
         exactScaleM012 = true,
         isBackground = true,
      },
      init = Plasma.Gyrokinetic.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z = xn[1], xn[2]
            local perturb = 1e-8*rho_e/L_T*math.cos(ky_min*y)
            return n0*(1+perturb)
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_T)
         end,
         exactScaleM012 = true,
      },
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkUpar", "GkTemp"}, 
   },

   -- adiabatic ions
   adiabaticIon = Plasma.AdiabaticSpecies {
      charge = qi,
      mass = mi,
      temp = Ti0,
      -- initial conditions
      init = function (t, xn, self)
         return n0
      end,
      evolve = false, -- evolve species?
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve fields?
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R0/(R0 + x)
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
print("expected growth rate = 2.92349*omega_de = ", 2.92349*omegade)
