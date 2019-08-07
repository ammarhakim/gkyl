-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- Maxwellian in 1x1v.
local function maxwellian1D(n, vx, ux, mass, temp)
   local v2 = (vx - ux)^2
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-mass*v2/(2*temp))
end

-- Normalization parameters.
epsilon0   = 1.0                       -- Permittivity of free space.
mu0        = 1.0                       -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- Speed of light.

elcMass   =  1.0         -- Electron mass.
elcCharge = -1.0         -- Electron charge.
ionMass   =  1836.153    -- Ion mass.
ionCharge =  1.0         -- Ion charge.

nuElc = 0.001                       -- Electron collision frequency.
nuIon = nuElc/math.sqrt(ionMass)    -- Ion collision frequency.

Te_Ti = 1.0    -- Ratio of electron to ion temperature.

n1  = 1.0
Te1 = 1.0
Ti1 = 1.0/Te_Ti

-- Derived parameters.
vtElc1 = math.sqrt(Te1/elcMass)
vtIon1 = math.sqrt(Ti1/ionMass)
-- Plasma frequency and Debye length.
wpe1 = math.sqrt(elcCharge^2*n1/(epsilon0*elcMass))
wpi1 = math.sqrt(ionCharge^2*n1/(epsilon0*ionMass))
lambdaD1 = vtElc1/wpe1

-- Domain size and simulation time.
LX = 64.0*lambdaD1

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 200.0/wpe1,       -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {0.0},            -- Configuration space lower left.
   upper       = {LX},             -- Configuration space upper right.
   cells       = {8},              -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   -- Boundary conditions for configuration space.
   periodicDirs = {1},    -- Periodic directions
   -- Decomposition for configuration space.
   decompCuts = {1},      -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.0025,

   -- Electrons.
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- Velocity space grid.
      lower = {-5.0*vtElc1},
      upper = {7.0*vtElc1},
      cells = {12},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vd = vtElc1
	 local xm = 0.25*LX
	 local fv = (1+math.exp(-0.75*(x-xm)^2))*maxwellian1D(n1, v, vd, elcMass, Te1) 
         if x>xm then
            fv = (1+math.exp(-0.075*(x-xm)^2))*maxwellian1D(n1, v, vd, elcMass, Te1)
         end
	 return fv
      end,
      evolve = true,    -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments           = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },

   -- Protons.
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- Velocity space grid
      lower = {-6.0*vtIon1+vtElc1},
      upper = { 6.0*vtIon1+vtElc1},
      cells = {12},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vd   = vtElc1
	 local xm   = 0.25*LX
	 local fv   = (1+math.exp(-0.75*(x-xm)^2))*maxwellian1D(n1, v, vd, ionMass, Ti1) 
         if x>xm then
            fv = (1+math.exp(-0.075*(x-xm)^2))*maxwellian1D(n1, v, vd, ionMass, Ti1)
         end
	 return fv
      end,
      evolve = true,    -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments           = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,    -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
