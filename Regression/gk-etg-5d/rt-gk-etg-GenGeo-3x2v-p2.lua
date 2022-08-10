-- Gkyl ------------------------------------------------------------------------
--
-- 5D GK ETG linear instability calculation with finite k_parallel and finite density gradient. 
-- parameters taken from Fig 2.4 of Beer's thesis.
--
--------------------------------------------------------------------------------
local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local math      = require("sci.math").generic

-- Physical parameters.
eV      = Constants.ELEMENTARY_CHARGE
qe      = -eV
qi      = eV
me      = Constants.ELECTRON_MASS
mi      = 2*Constants.PROTON_MASS   -- Deuterium ions.
Te0     = 2072*eV 
Ti0     = 2072*eV 
B0      = 1.91   -- maximum field strength on axis [T]
R0      = 1.313  -- major radius of magnetic axis [m]
a       = 0.4701 -- minor radius [m]
n0      = 4.992*10^(19) -- [1/m^3]
eps_n   = 0.2
eta_e   = 3.0
ky_rhoe = 1.2
kz_Ln   = 0.1
-- Derived parameters.
r0         = 0.5*a     -- Minor radius of flux tube.
R          = R0 + r0   -- Major radius of flux tube.
B          = B0*R0/R   -- Magnetic field strength in flux tube.
vte  	   = math.sqrt(Te0/me)
c_s        = math.sqrt(Te0/mi)
omega_ci   = math.abs(qi*B/mi)
omega_ce   = math.abs(qe*B/me)
rho_s      = c_s/omega_ci
rho_e      = vte/omega_ce
ky_min     = ky_rhoe / rho_e
dr         = 2*math.pi/ky_min
L_n        = R*eps_n
L_T        = L_n/eta_e
kz_min     = kz_Ln / L_n
L_parallel = 2*math.pi/kz_min
-- Velocity grid parameters.
N_VPAR, N_MU = 16, 8
VPAR_UPPER = math.min(4, 2.5*math.sqrt(N_VPAR/4))*vte
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = math.min(16, 8*math.sqrt(N_MU/2))*me*vte*vte/B/2
omegade  = ky_min*rho_e*vte/R

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = .15e-8, -- End time.
   nFrame = 1,       -- Number of output frames.
   lower  = {r0 - dr/2, -dr/2, -L_parallel/2}, -- Configuration space lower left.
   upper  = {r0 + dr/2,  dr/2, L_parallel/2},  -- Configuration space upper right.
   mapc2p = function(xc)
      -- Field-aligned coordinates (x,y,z).
      local x, y, z = xc[1], xc[2], xc[3]
      -- Cylindrical coordinates (R,phi).
      local R = x+R0
      local phi = z/(R0+r0)
      -- Cartesian coordinates (X,Y,Z).
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      local Z = y
      return X, Y, Z
   end,
   cells       = {2, 6, 6},     -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2" or "rk3".
   cflFrac     = 1.0,

   -- Decomposition for configuration space.
   decompCuts = {1, 1, 1}, -- Cuts in each configuration direction.

   -- Poundary conditions for configuration space.
   periodicDirs = {1,2,3}, -- Periodic directions.

   deltaF = true,

   -- Gyrokinetic electrons.
   electron = Plasma.Species {
      charge = qe,
      mass = me,
      -- Velocity space grid.
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      decompCuts = {1, 1},
      -- Initial conditions.
      background = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x = xn[1]
            return n0*(1-(x-r0)/L_n)
         end,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_T)
         end,
      },
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            local perturb = 1e-5*rho_e/L_T*math.cos(ky_min*y+kz_min*z)
            return n0*(1-(x-r0)/L_n) + n0*perturb
         end,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_T)
         end,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "M2"},
   },

   -- Adiabatic ions
   adiabaticIon = Plasma.AdiabaticSpecies {
      charge = qi,
      mass = mi,
      temp = Ti0,
      -- Initial conditions.
      init = function (t, xn)
         local x = xn[1]
         return n0*(1-(x-r0)/L_n)
      end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
      --polarizationWeight = 0.0,
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return B0*R0/(R0 + x)
      end,
      geometryType = "GenGeo",
      evolve       = false,   -- Geometry is not time-dependent.
   },
}
-- Run application.
plasmaApp:run()

print("expected growth rate = ", 2.26579*vte/R)
