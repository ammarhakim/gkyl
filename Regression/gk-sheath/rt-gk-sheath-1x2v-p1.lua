-- Gkyl ------------------------------------------------------------------------
--
-- This is a 1x2v test based off a 3x2v NSTX-like SOL simulation.
--
--------------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"

-- Universal constant parameters.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
qe, qi   = -eV, eV
me, mp   = Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters.
mi  = 2.014*mp   -- Deuterium ions.
Te0 = 40*eV 
Ti0 = 40*eV 
n0  = 7e18       -- [1/m^3]

-- Geometry and magnetic field.
B_axis = 0.5    -- [T]
R0     = 0.85   -- [m]
a0     = 0.15   -- [m]
R      = R0 + a0
B0     = B_axis*(R0/R) -- [T]
Lpol   = 2.4    -- [m]

-- Perpendicular wavenumber (k_perp) times ion sound gyroradius (rho_s).
kperpRhos = 0.3

-- Parameters for collisions.
nuFrac = 0.1
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))

-- Derived parameters
vti, vte = math.sqrt(Ti0/mi), math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci

-- Box size.
Lx = 50*rho_s
Ly = 100*rho_s
Lz = 4 -- [m]

-- Source parameters.
P_SOL   = 3.4e6 -- [W], total SOL power, from experimental heating power
P_src   = P_SOL*Ly*Lz/(2*math.pi*R*Lpol) -- [W], fraction of total SOL power into flux tube
xSource = R     -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature

-- Source profiles.
sourceDensity = function (t, xn)
   local z = xn[1]
   local sourceFloor = 0.1
   if math.abs(z) < Lz/4 then
      return math.max(1., sourceFloor)
   else
      return 1e-40
   end
end
sourceTemperature = function (t, xn)
   local z = xn[1]
   return 80*eV
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 6.e-6,                   -- End time.
   nFrame      = 1,                       -- Number of output frames.
   lower       = {-Lz/2},                  -- Configuration space lower left.
   upper       = { Lz/2},                  -- Configuration space upper right.
   cells       = {8},                      -- Configuration space cells.
   basis       = "serendipity",            -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                        -- Polynomial order.
   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.4,
   restartFrameEvery = .5,
   calcIntQuantEvery = 1./60.,

   decompCuts  = {1},    -- MPI subdomains/processes.

   -- Gyrokinetic electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      lower = {-4*vte, 0},
      upper = { 4*vte, 12*me*vte^2/(2*B0)},
      cells = {6, 4},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local z = xn[1]
            local Ls              = Lz/4
            local effectiveSource = sourceDensity(t,{0})
            local c_ss            = math.sqrt(5/3*sourceTemperature(t,{0})/mi)
            local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
            local perturb         = 0 
            if math.abs(z) <= Ls then
               return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
            else
               return nPeak/2*(1+perturb)
            end
         end,
         temperature = function (t, xn)
            local z = xn[1]
            return Te0
         end,
         scaleWithSourcePower = true,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
      source = Plasma.Source {
         kind        = "Maxwellian",
         density     = sourceDensity,
         temperature = sourceTemperature,
         power       = P_src/2,
         diagnostics = {"M0","Temp","intM0","intM2"},
      },
      polarizationDensityFactor = n0,
      evolve = true, -- Evolve species?
      --applyPositivity = true,
      diagnostics = {"M0", "Upar", "Temp", "Beta", "intM0", "intM1", "intEnergy",}, 
      randomseed = randomseed,
      bcx = {Plasma.SheathBC{diagnostics={"M0","Upar","Temp","Beta","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.SheathBC{diagnostics={"M0","Upar","Temp","Beta","Energy","intM0","intM1","intKE","intEnergy"}}},
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      charge = qi,  mass = mi,
      -- Velocity space grid.
      lower = {-4*vti, 0},
      upper = { 4*vti, 12*mi*vti^2/(2*B0)},
      cells = {6, 4},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local z = xn[1]
            local Ls              = Lz/4
            local effectiveSource = sourceDensity(t,{0})
            local c_ss            = math.sqrt(5/3*sourceTemperature(t,{0})/mi)
            local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
            local perturb         = 0
            if math.abs(z) <= Ls then
               return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
            else
               return nPeak/2*(1+perturb)
            end
         end,
         temperature = function (t, xn)
            local z = xn[1]
            return Ti0
         end,
         scaleWithSourcePower = true,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
      source = Plasma.Source {
         kind        = "Maxwellian",
         density     = sourceDensity,
         temperature = sourceTemperature,
         power       = P_src/2,
         diagnostics = {"M0","Temp","intM0","intM2"},
      },
      polarizationDensityFactor = n0,
      evolve = true, -- Evolve species?
      --applyPositivity = true,
      diagnostics = {"M0", "Upar", "Temp", "intM0", "intM1", "intKE", "intEnergy"}, 
      randomseed = randomseed,
      bcx = {Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0","intM1","intKE","intEnergy"}}},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
      kperpSq = (kperpRhos/rho_s)^2,
      -- Dirichlet in x, periodic in y. Potential phi has homogeneous Neumann
      -- BC for the smoothing operation that enforces continuity in z.
      bcLowerPhi = {{T = "N", V = 0.0}},
      bcUpperPhi = {{T = "N", V = 0.0}},
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local z = xn[1]
         return B0
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
