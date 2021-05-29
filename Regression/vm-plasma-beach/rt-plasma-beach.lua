-- Gkyl ------------------------------------------------------------------------
local Plasma    = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Constants = require "Lib.Constants"

-- Constants.
epsilon0   = Constants.EPSILON0
mu0        = Constants.MU0
lightSpeed = Constants.SPEED_OF_LIGHT

eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2.014*Constants.PROTON_MASS   -- Deuterium ions.

Te_Ti   = 1.0
elcTemp = 1.0*Constants.EV2KELVIN
ionTemp = elcTemp/Te_Ti

-- Thermal speeds.
vtElc = math.sqrt(Constants.BOLTZMANN_CONSTANT*elcTemp/me)
vtIon = math.sqrt(Constants.BOLTZMANN_CONSTANT*ionTemp/mi)

-- Domain size and simulation time.
xlower = 0.0   -- Lower bounds of domain.
xupper = 1.0   -- Upper bounds of domain.

LX = xupper-xlower

-- Resolution, time-stepping etc.
NX = 24
NVX = 6
NVY = 6

tEnd = 5.0e-9
nFrames = 1

-- Compute coordinate of interior last edge.
dx        = (xupper-xlower)/NX
xLastEdge = xupper-dx

dx100 = (xupper-xlower)/100
-- Compute drive frequency.
deltaT     = dx100/lightSpeed
driveOmega = math.pi/10/deltaT

-- Maxwellian in 1x2v.
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = tEnd,            -- End time.
   nFrame      = nFrames,         -- Number of output frames.
   lower       = { 0.0 },         -- Configuration space lower left
   upper       = { LX },          -- Configuration space upper right.
   cells       = {NX},            -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2" or "rk3".
   cflFrac     = 1/25.0,
   -- Decomposition for configuration space.
   decompCuts = {1},     -- Cuts in each configuration direction.
   useShared  = false,   -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {},   -- Periodic directions.

   -- Electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      -- Velocity space grid
      lower = {-6.0*vtElc, -6.0*vtElc},
      upper = { 6.0*vtElc,  6.0*vtElc},
      cells = {NVX, NVY},
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         local wpdt   = 25*(1-x)^5   -- Plasma frequency
         local factor = deltaT^2*qi*qi/(me*epsilon0)
         local nElc   = wpdt^2/factor

	 return maxwellian2D(nElc, vx, vy, 0.0, 0.0, vtElc)
      end,
      -- Boundary conditions.
      bcx    = {Plasma.Species.bcCopy, Plasma.Species.bcCopy},
      evolve = true,   -- Evolve species?

      diagnostics = { "M0", "M1i", "M2" }
   },

   -- Protons.
   ion = Plasma.Species {
      charge = qi,  mass = mi,
      -- Velocity space grid.
      lower = {-6.0*vtIon, -6.0*vtIon},
      upper = { 6.0*vtIon,  6.0*vtIon},
      cells = {NVX, NVY},
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]

         local wpdt   = 25*(1-x)^5   -- Plasma frequency
         local factor = deltaT^2*qi*qi/(me*epsilon0)
         local nIon   = wpdt^2/factor

	 return maxwellian2D(nIon, vx, vy, 0.0, 0.0, vtIon)
      end,
      -- Boundary conditions.
      bcx    = {Plasma.Species.bcCopy, Plasma.Species.bcCopy},
      evolve = false,   -- Evolve species?

      diagnostics = { "M0", "M1i", "M2" }
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
	 local x = xn[1]
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      -- Boundary conditions.
      bcx    = {Plasma.Field.bcCopy, Plasma.Field.bcCopy},
      evolve = true,   -- Evolve field?
   },

   -- Current antenna.
   driveSpecies = Plasma.FuncSpecies {
      charge = 1.0,  mass = 1.0,
      momentumDensity = function (t, xn)
	 local x  = xn[1]
         local J0 = 1.0   -- Amps/m^3.
         local Jy = 0.0
         if (x>xLastEdge) then Jy = -J0*math.sin(driveOmega*t) end
	 return 0.0, Jy
      end,
      evolve = true,   -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
