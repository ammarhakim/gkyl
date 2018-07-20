-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- physical parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/100.0
lightSpeed = 1.0
epsilon0 = 1.0
mu0 = 1.0
mgnErrorSpeedFactor = 1.0

n0 = 1.0
VAe = 0.5
plasmaBeta = 1.0
lambda = 0.5
TiOverTe = 5.0
nbOverN0 = 0.2
pert = 0.1
Valf = VAe*math.sqrt(elcMass/ionMass)
B0 = Valf*math.sqrt(n0*ionMass)
OmegaCi0 = ionCharge*B0/ionMass
psi0 = pert*B0

-- domain size
Lx = 25.6
Ly = 12.8

momentApp = Plasma.App {
   logToFile = true,

   tEnd = 40.0,
   nFrame = 1,
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {64, 32},
   timeStepper = "fvSplit",

   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Plasma.Moments.Species {
      charge = -1.0, mass = 1.0,
      -- initial conditions
      init = function (t, xn)
	 return 1.0, 0.0, 0.0, 0.0, 1.0
      end,
      evolve = true, -- evolve species?

   },
}
-- run application
momentApp:run()
