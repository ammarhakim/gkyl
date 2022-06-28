-- Gkyl ------------------------------------------------------------------------
--
-- 4D GK ETG linear instability calculation
-- using 'pgkyl -f etg4d_elecEnergy_ growth' should approximately give growth rate printed at end of run
--
--------------------------------------------------------------------------------
local Plasma    = require ("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local math      = require("sci.math").generic

-- Physical parameters.
eV  = Constants.ELEMENTARY_CHARGE
qe  = -eV
qi  = eV
me  = Constants.ELECTRON_MASS
mi  = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0 = 2072*eV 
Ti0 = 2072*eV 
B0  = 1.91     -- Magnetic field strength on axis [T].
R0  = 1.313    -- Major radius of magnetic axis [m].
a   = 0.4701   -- Minor radius [m].
n0  = 4.992e19 -- [1/m^3].
-- Derived parameters.
r0       = 0.5*a     -- Minor radius of flux tube.
R        = R0 + r0   -- Major radius of flux tube.
B        = B0*R0/R   -- Magnetic field strength in flux tube.
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
-- Velocity grid parameters.
N_VPAR, N_MU = 16, 8
VPAR_UPPER   = 4*vte
VPAR_LOWER   = -VPAR_UPPER
MU_LOWER     = 0
MU_UPPER     = 16*me*vte*vte/B/2

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = .5e-6, -- End time.
   nFrame = 2,     -- Number of output frames.
   lower  = {r0 - 0.001*dr/2, -dr/2}, -- Configuration space lower left.
   upper  = {r0 + 0.001*dr/2,  dr/2}, -- Configuration space upper right.
   cells  = {1, 8},   -- Configuration space cells.
   mapc2p = function(xc)
      local r, y = xc[1], xc[2]
      local X = r+R0
      local Y = y
      return X, Y
   end,
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2" or "rk3".
   cflFrac     = .9,

   -- Boundary conditions for configuration space.
   periodicDirs = {1,2}, -- Periodic directions.
   deltaF       = true,  -- Only apply BCs to fluctuations, and use perturbed moments in field solve.

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
            return n0
         end,
         driftSpeed = 0.0,
         temperature = function (t, xn)
            local x = xn[1]
            return Te0*(1-(x-r0)/L_T)
         end,
         exactScaleM012 = true,
      },
      init = Plasma.MaxwellianProjection {
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
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "intM0", "intM2"}, 
   },

   -- Adiabatic ions.
   adiabaticIon = Plasma.AdiabaticSpecies {
      charge = qi,
      mass = mi,
      temp = Ti0,
      -- Initial conditions.
      init = function (t, xn, self)
         return n0
      end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R0/(R0 + x)
      end,
      evolve = false,   -- Geometry is not time-dependent.
   },
}
-- Run application.
plasmaApp:run()
--print("expected growth rate = 2.92349*omega_de = ", 2.92349*omegade)
