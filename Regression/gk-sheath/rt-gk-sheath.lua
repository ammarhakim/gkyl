-- This test is based off NSTX-like SOL simulation
-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"
local Mpi = require "Comm.Mpi"

-- Universal constant parameters.
eps0         = Constants.EPSILON0
eV           = Constants.ELEMENTARY_CHARGE
qe           = -eV
qi           = eV
me           = Constants.ELECTRON_MASS

-- Plasma parameters.
mi           = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0          = 40*eV 
Ti0          = 40*eV 
n0           = 7e18  -- [1/m^3]

-- Geometry and magnetic field.
B_axis       = 0.5   -- [T]
R0           = 0.85  -- [m]
a0           = 0.15   -- [m]
R            = R0 + a0
B0           = B_axis*(R0/R) -- [T]

-- Source parameters.
P_SOL        = 8.1e5 -- [W] 
S0           = 5.7691e23
xSource      = R -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature

-- Parameters for collisions.
nuFrac = 1.0
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))

-- Derived parameters
vti      = math.sqrt(Ti0/mi)
vte  	 = math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci

-- Box size.
Lx = 50*rho_s
Ly = 100*rho_s
Lz = 4 -- [m]

-- Source profiles.
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local sourceFloor = 0.1
   if math.abs(z) < Lz/4 then
      return 0.90625*S0*math.max(math.exp(-(x-xSource)^2/(2*lambdaSource)^2), sourceFloor)
   else
      return 1e-10
   end
end
sourceTemperature = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < xSource + 3*lambdaSource then
      return 80*eV
   else
      return 30*eV
   end
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 1e-6,                     -- End time.
   nFrame      = 1,                     -- Number of output frames.
   lower       = {R - Lx/2, -Ly/2, -Lz/2}, -- Configuration space lower left.
   upper       = {R + Lx/2, Ly/2, Lz/2},   -- Configuration space upper right.
   cells       = {4, 1, 8},              -- Configuration space cells.
   basis       = "serendipity",            -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                        -- Polynomial order.
   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.3,
   restartFrameEvery = .1,

   -- Decomposition for configuration space.
   decompCuts = {1, 1, 1}, -- Cuts in each configuration direction.
   useShared = false,      -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {2},     -- Periodic in y only.

   -- Gyrokinetic electrons.
   electron = Plasma.GkSpecies {
      charge = qe,
      mass  = me,
      lower = {-4*vte, 0},
      upper = {4*vte, 12*me*vte^2/(2*B0)},
      cells = {8, 4},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = {"maxwellian", 
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local Ls              = Lz/4
                 local effectiveSource = sourceDensity(t,{x,y,0})
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb         = 0 
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if (x < xSource + 3*lambdaSource) then 
                    return 50*eV
                 else 
                    return 20*eV
                 end
              end,
      },
      --coll   = Plasma.GkLBOCollisions { collFreq = nuElc },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperature},
      evolve = true, -- Evolve species?
      --applyPositivity = true,
      diagnosticMoments = {"GkM0", "GkUpar", "GkTemp"}, 
      randomseed = randomseed,
      bcx = {Plasma.GkSpecies.bcZeroFlux, Plasma.GkSpecies.bcZeroFlux},
      bcz = {Plasma.GkSpecies.bcSheath, Plasma.GkSpecies.bcSheath},
   },

   -- Gyrokinetic ions
   ion = Plasma.GkSpecies {
      charge = qi,
      mass   = mi,
      -- Velocity space grid.
      lower = {-4*vti, 0},
      upper = {4*vti, 12*mi*vti^2/(2*B0)},
      cells = {8, 4},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = {"maxwellian", 
              density = function (t, xn)
                 local x, y, z         = xn[1], xn[2], xn[3]
                 local Ls              = Lz/4
                 local effectiveSource = sourceDensity(t,{x,y,0})
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb         = 0
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if x < xSource + 3*lambdaSource then 
                    return 50*eV
                 else 
                    return 20*eV
                 end
              end,
      },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperature},
      evolve = true, -- Evolve species?
      diagnosticMoments = {"GkM0", "GkUpar", "GkTemp"}, 
      randomseed = randomseed,
      bcx = {Plasma.GkSpecies.bcZeroFlux, Plasma.GkSpecies.bcZeroFlux},
      bcz = {Plasma.GkSpecies.bcSheath, Plasma.GkSpecies.bcSheath},
   },

   -- Field solver.
   field = Plasma.GkField {
      -- Dirichlet in x.
      phiBcLeft  = { T ="D", V = 0.0},
      phiBcRight = { T ="D", V = 0.0},
      -- Periodic in y. --
      -- No bc in z.
      phiBcBack  = { T ="N", V = 0.0},
      phiBcFront = { T ="N", V = 0.0},
      evolve     = true, -- Evolve fields?
   },

   -- Magnetic geometry.
   funcField = Plasma.GkGeometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R/x
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- run application
plasmaApp:run()
