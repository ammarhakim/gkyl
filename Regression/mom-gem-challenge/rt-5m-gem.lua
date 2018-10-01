-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

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

momentApp = Moments.App {
   logToFile = true,

   tEnd = 40.0,
   nFrame = 1,
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {64, 32},
   timeStepper = "fvDimSplit",

   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local TeFrac = 1.0/(1.0 + TiOverTe)
	 local sech2 = (1.0/math.cosh(y/lambda))^2
	 local n = n0*(sech2 + nbOverN0)
	 local Jz = -(B0/lambda)*sech2
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoe = n*elcMass
	 local ezmom = (elcMass/elcCharge)*Jz*TeFrac
	 local ere = n*Ttotal*TeFrac/(gasGamma-1) + 0.5*ezmom*ezmom/rhoe
	 
	 return rhoe, 0.0, 0.0, ezmom, ere
      end,
      evolve = true, -- evolve species?
      
   -- bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   bcy = { Moments.Species.bcReflect, Moments.Species.bcReflect },
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local TiFrac = TiOverTe/(1.0 + TiOverTe)
	 local sech2 = (1.0/math.cosh(y/lambda))^2
	 local n = n0*(sech2 + nbOverN0)
	 local Jz = -(B0/lambda)*sech2
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoi = n*ionMass
	 local izmom = (ionMass/ionCharge)*Jz*TiFrac
	 local eri = n*Ttotal*TiFrac/(gasGamma-1) + 0.5*izmom*izmom/rhoi
	 
	 return rhoi, 0.0, 0.0, izmom, eri
      end,
      evolve = true, -- evolve species?
      
   -- bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   bcy = { Moments.Species.bcReflect, Moments.Species.bcReflect },
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local alpha = perturbation
	 local k = knumber
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   -- bcx = { Moments.Field.bcCopy, Moments.Field.bcCopy },
   bcy = { Moments.Field.bcReflect, Moments.Field.bcReflect },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "time-centered",
   },   

}
-- run application
momentApp:run()
