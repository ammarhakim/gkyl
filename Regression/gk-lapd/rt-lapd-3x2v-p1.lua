-- Gkyl ------------------------------------------------------------------------
--
-- Reproduce Eric Shi's LAPD case.
--
-- E. Shi, et al. JPP 83, 905830304 (2017).
--
--------------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"
local xsys      = require "xsys"
local Logger    = require "Lib.Logger"

local log = Logger {
   logToFile = xsys.pickBool(logToFile, true)
}

-- Universal constants.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
qe, qi   = -eV, eV

mRat = 400
mi   = 3.973*Constants.PROTON_MASS -- He ion mass (kg).
me   = mi/mRat -- Electron mass (kg).
Te0  = 6*eV    -- Electron temperature (J).
Ti0  = 1*eV    -- Ion temperature (J).
B0   = 0.0398  -- Magnetic field magnitude (T).
n0   = 2.e18   -- Particle density (1/m^3).

-- Derived parameters.
vti      = math.sqrt(Ti0/mi)
vte      = math.sqrt(Te0/me) -- g1 uses factor of .5 in sqrt arg
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci

-- Simulation box size (m).
Lx = 100.*rho_s
Ly = 100.*rho_s
R  = 40.*rho_s
Lz = 36.*R       -- approx 18.0492 m.

-- Grid parameters. 
Nz = 8          -- Number of cells in z.
pOrder = 1       -- Polynomial order of DG basis.

-- Bias parameters.
V_bias = 0   -- Biasing voltage (V).
r_bias = .26 -- Radius of the biased plate (m).

-- Source parameters.
L_s = 0.5*rho_s
r_s = 20*rho_s
S0  = 1.08*n0*c_s/Lz

-- Collision parameters.
nuFrac = 0.1
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc     = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))  --collision freq

logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon     = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))

nuElcIon = math.sqrt(2)*nuElc
nuIonElc = nuElcIon/mRat

log(string.format("Reference parameters (SI units):\n"))
log(string.format("  mi     = %8.4e\n",mi ))
log(string.format("  me     = %8.4e\n",me ))
log(string.format("  Te0    = %8.4e\n",Te0))
log(string.format("  Ti0    = %8.4e\n",Ti0))
log(string.format("  n0     = %8.4e\n",n0 ))
log(string.format("  B0     = %8.4e\n",B0 ))
log(string.format("  rho_s  = %8.4e\n",rho_s))
log(string.format("  V_bias = %8.4e\n",V_bias))
log(string.format("  r_bias = %8.4e\n",r_bias))
log(string.format("  r_s    = %8.4e\n",r_s   ))
log(string.format("  L_s    = %8.4e\n",L_s   ))
log(string.format("  S0     = %8.4e\n",S0    ))
log(string.format("  Lx     = %8.4e\n",Lx))
log(string.format("  Ly     = %8.4e\n",Ly))
log(string.format("  Lz     = %8.4e\n",Lz))
log(string.format("\n"))

-- Compute the omega_H mode frequency. The time step has to stay below it.
local dzEff    = Lz/(Nz*(pOrder+1))
local kParMax  = 2.*math.pi/(2.*dzEff)
local kperpMin = 2.*math.pi/(2.*Lx)
local omega_H  = kParMax*vte/(kperpMin*rho_s)

log(string.format("  Max omega_H = kpar_max*vte/(kperp_min*rhos) = %e rad/s\n", omega_H))
log("\n")

-- Initial profiles.
local profileA = function(r, c_edge)
   local Lperp = Lx
   if r < Lperp/2. then
      return (1.-c_edge)*(1.-(r/(Lperp/2.))^2)^3 + c_edge
   else
      return c_edge
   end
end
local initDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local r       = math.sqrt(x^2 + y^2)
   return n0*profileA(r, 1./20.)
end
local initTemperatureElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local r       = math.sqrt(x^2 + y^2)
   local TeInit  = 5.7*eV
   return TeInit*profileA(r, 1./5.)
end
local initTemperatureIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 1.*eV
end

-- Source profiles.
local sourceDensity = function (t, xn)
   local x, y, z     = xn[1], xn[2], xn[3]
   local r           = math.sqrt(x^2 + y^2)
   local sourceFloor = 0.01
   return S0*(sourceFloor + (1.-sourceFloor)*0.5*(1 - math.tanh((r-r_s)/L_s)))
end
local sourceTemperatureElc = function (t, xn)                                                                       
   local x, y, z = xn[1], xn[2], xn[3]
   local r       = math.sqrt(x^2 + y^2)
   local TeS     = 6.8*eV
   return TeS*profileA(r, 1./2.5)
end
local sourceTemperatureIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 1.*eV
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD) --+os.time()

local finalTime = .25e-6
local numFrames = 1

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = finalTime,             -- End time.
   nFrame      = numFrames,             -- Number of output frames.
   lower       = {-Lx/2, -Ly/2, -Lz/2}, -- Configuration space lower left.
   upper       = { Lx/2,  Ly/2,  Lz/2}, -- Configuration space upper right.
   cells       = {18, 18, Nz},          -- Configuration space cells.
   basis       = "serendipity",         -- One of "serendipity" or "maximal-order".
   polyOrder   = pOrder,                -- Polynomial order.
   timeStepper = "rk3",                 -- One of "rk2" or "rk3".
--   cflFrac     = 0.2,
   restartFrameEvery = 0.05,
   calcIntQuantEvery = 1./(10.*numFrames),  -- Aim for 10x more frequently than frames.
   groupDiagnostics  = true,

   -- Decomposition for configuration space.
   decompCuts = {1, 1, 1},             -- Cuts in each configuration direction.
   parallelizeSpecies = true,

   -- Gyrokinetic electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      -- Velocity space grid.
      lower = {-4.*vte, 0.},
      upper = { 4.*vte, 0.75*me*(4.*vte)^2/(2.*B0)},
      cells = {10, 5},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            local perturb = 1e-3*(1.-0.5*math.random())
	    return initDensity(t,xn)*(1.+perturb)
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return initTemperatureElc(t,xn)
         end,
      },
      coll = Plasma.LBOCollisions {
	 collideWith = {'elc', 'ion'},
         frequencies = {nuElc, nuElcIon},
--         normNu      = {nuElc*((2*(vte^2))^(3/2))/n0, nuElcIon*((vte^2+vti^2)^(3/2))/n0},
      },
      polarizationDensityFactor = n0,
      source     = Plasma.Source{density = sourceDensity, temperature = sourceTemperatureElc},
      randomseed = randomseed,
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcy = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{}, Plasma.SheathBC{}},
      evolve = true,   -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "intM0"}, 
      nDistFuncFrame = 10, 
   },

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = qi,  mass = mi,
      -- Velocity space grid.
      lower = {-4.*vti, 0},
      upper = { 4.*vti, 0.75*mi*(4.*vti)^2/(2.*B0)},
      cells = {10, 5},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            local perturb = 1e-3*(1.-0.5*math.random())
	    return initDensity(t,xn)*(1.+perturb)
         end,
         temperature = function (t, xn)
            local x = xn[1]
	    return initTemperatureIon(t,xn)
         end,
      },
      coll = Plasma.LBOCollisions {
	 collideWith = {'ion', 'elc'},
         frequencies = {nuIon, nuIonElc},
--         normNu      = {nuIon*((2*(vti^2))^(3/2))/n0, nuIonElc*((vti^2+vte^2)^(3/2))/n0},
      },
      polarizationDensityFactor = n0,
      source     = Plasma.Source{density = sourceDensity, temperature = sourceTemperatureIon},
      randomseed = randomseed,
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcy = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{}, Plasma.SheathBC{}},
      evolve = true,   -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "intM0"}, 
      nDistFuncFrame = 10, 
   },

   -- Field solver.
   field = Plasma.Field {
      -- Dirichlet in x, Dirichlet in y, no bc in z.
      bcLowerPhi = {{T ="D", V = 0.0}, {T ="D", V = 0.0}, {T ="N", V = 0.0}},
      bcUpperPhi = {{T ="D", V = 0.0}, {T ="D", V = 0.0}, {T ="N", V = 0.0}},
      evolve = true,   -- Evolve fields?
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         return B0
      end,

      phiWall = function(t, xn)
         local x, y = xn[1], xn[2]
	 local r    = math.sqrt(x^2 + y^2)
         if r < r_s then
            return V_bias
         else
            return 0
         end
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()

