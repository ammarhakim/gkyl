-- Gkyl ------------------------------------------------------------------------
--
-- An 3x2v simulation of the LTX edge including open and closed field lines,
-- a artificial q profile that matches TRANSP at the LCFS, and uses 
-- parameters similar to those near the LCFS of LTX's shot 103795
-- at t=469.11 ms.
--
-- The plasma parameters are taken from Anurag Maan & Elizabeth Perez's fresh
-- Li simulation.
--
--------------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"
local xsys      = require "xsys"
local Logger    = require "Lib.Logger"
local math      = require("sci.math").generic
local quad      = require("sci.quad")
local diff      = require("sci.diff-recursive")
local df        = diff.df

local log = Logger { logToFile = xsys.pickBool(logToFile, true) }

-- Universal constant parameters.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
qe, qi   = -eV, eV
me, mp   = Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters. Chosen based on the value of a cubic sline
-- between the last TS data inside the LCFS and the probe data in
-- in the far SOL, near R=0.475 m.
mi  = mp   -- Hydrogen ions.
Te0 = 250*eV 
Ti0 = 250*eV 
n0  = 3.0e18     -- [1/m^3]

-- Geometry and magnetic field.
R_axis     = 0.406051632        -- [m]
B_axis     = 0.240606108        -- [T]
R_LCFSmid  = 0.61183           -- Major radius of the LCFS at the outboard midplane [m].
Rmid_min   = R_LCFSmid-2*0.05457   -- Minimum midplane major radius of simulation box [m].
Rmid_max   = R_LCFSmid+0.05457   -- Maximum midplane major radius of simulation box [m].
R0         = 0.5*(Rmid_min+Rmid_max)    -- Major radius of the simulation box [m].
a_mid      = R_LCFSmid-R_axis   -- Minor radius at outboard midplane [m].
r0         = R0-R_axis          -- Minor radius of the simulation box [m].
B0         = B_axis*(R_axis/R0) -- Magnetic field magnitude in the simulation box [T].
q_LCFS     = 4.3153848          -- Safety factor at the LCFS.
s_LCFS     = 2.6899871          -- Magnetic shear at the LCFS. Should be ~6.6899871 but that makes qprofile change sign in the SOL.
kappa      = 1.3                -- Elongation (=1 for no elongation).
delta      = 0.4                -- Triangularity (=0 for no triangularity).

function r_x(x) return Rmid_min+x-R_axis end   -- Minor radius given x, where x \in [0,Lx].
function qprofile(r) return q_LCFS/(1.-s_LCFS*((r-a_mid)/a_mid)) end   -- Magnetic safety profile.

q0 = qprofile(r0)    -- Magnetic safety factor in the center of domain.

log(string.format(" Rmid_min = %f\n",Rmid_min))
log(string.format(" Rmid_max = %f\n",Rmid_max))
log(string.format(" a_mid    = %f\n",a_mid   ))
log(string.format(" R0       = %f\n",R0      ))
log(string.format(" r0       = %f\n",r0      ))
log(string.format(" q0       = %f\n",q0      ))
log(string.format(" B0       = %f\n",B0      ))

nuFrac = 1.0
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))
-- Electron-ion and ion-electron collision frequencies.
nuElcIon = nuElc*math.sqrt(2)
nuIonElc = nuElcIon/(mi/me)

-- Derived parameters
vti, vte = math.sqrt(Ti0/mi), math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci
log(string.format(" rho_s = %f\n",rho_s))

-- Box size.
Lx = Rmid_max-Rmid_min
xMin, xMax = 0., Lx
rMin, rMax = Rmid_min-R_axis, Rmid_max-R_axis
ntoroidal = 3.0
Ly = 2.*math.pi*(r0/q0)/ntoroidal
Lz = 2.*math.pi
x_LCFS = R_LCFSmid - Rmid_min  -- Location of the LCFS in x.
log(string.format(" ntoroidal = %f\n",ntoroidal))
log(string.format(" Lx = %f\n",Lx))
log(string.format(" Ly = %f\n",Ly))
log(string.format(" Lz = %f\n",Lz))
log(string.format(" x_LCFS = %f\n",x_LCFS))

epsilon0  = r0/R0              -- Inverse aspect ratio in the center of the domain.
nuStarElc = nuElc*q0*R0/(vte*(epsilon0^(3./2.)))
nuStarIon = nuIon*q0*R0/(vti*(epsilon0^(3./2.)))
log(string.format(" nuStarElc = %g\n", nuStarElc))
log(string.format(" nuStarIon = %g\n", nuStarIon))


-- Functions needed for (x,y,z) -> (X,Y,Z) mapping in Miller geometry.
local function R(r, theta) return R_axis + r*math.cos(theta + math.asin(delta)*math.sin(theta)) end
local function Z(r, theta) return kappa*r*math.sin(theta) end
local function Bphi(R) return B0*R0/R end
local function Jr(r, theta)
   return R(r,theta)*(df(R,1)(r,theta)*df(Z,2)(r,theta) - df(Z,1)(r,theta)*df(R,2)(r,theta))
end
local function dPsidr(r, theta)
   local function integrand(t) return Jr(r,t)/R(r,t)^2 end
   local integral
   integral, _ = quad.dblexp(integrand, 0, 2*math.pi, 1e-10)
   return B0*R_axis/(2*math.pi*qprofile(r))*integral
end
local function J(r, theta) return Jr(r, theta)/dPsidr(r, theta) end
local function alpha(r, theta, phi)
   local function integrand(t) return Jr(r,t)/R(r,t)^2 end
   local integral
   local t = theta
   while diff.lt(t, -math.pi) do t = t+2*math.pi end
   while diff.lt( math.pi, t) do t = t-2*math.pi end
   if diff.lt(0, t) then
      integral, _ =  quad.dblexp(integrand, 0, t, 1e-10)/dPsidr(r,theta)
   else
      integral, _ = -quad.dblexp(integrand, t, 0, 1e-10)/dPsidr(r,theta)
   end
   return phi - B0*R_axis*integral
end
local function gradr(r, theta)
   return R(r,theta)/Jr(r,theta)*math.sqrt(df(R,2)(r,theta)^2 + df(Z,2)(r,theta)^2)
end

--START the canonical section
local function psi_fit(r)
   local a, b, c, d, e =-1.1463112730315916e-02,  3.2553352140788654e+00, -1.6917998161853689e+01, -1.4622473121061603e+00, 3.2671388745728050e-03
   return a/(b + math.exp(-c*r + d)) + e
end

local BtMax = Bphi(R(rMax,Lz/2))
local BpMax = dPsidr(rMax,Lz/2)/R(rMax,Lz/2)*gradr(rMax,Lz/2)
local BMax = math.sqrt(BtMax^2 + BpMax^2)

local function canonical_maxwellian(xn,q_s,m_s,T0)
   local x,y,z,vpar,mu = xn[1], xn[2], xn[3], xn[4], xn[5]
   local r = r_x(x)
   local psi0 = psi_fit(Rmid_min-R_axis)
   local psi = psi_fit(r)
   local Bt = Bphi(R(r,z))
   local Bp = dPsidr(r,z)/R(r,z)*gradr(r,z)
   local B = math.sqrt(Bt^2 + Bp^2)
   --Calculate the Density and temperature profiles with x -> xCanon
   --xCanon(x) is given by \bar{r}(r) in Dif-Pradalier. Note phi = 0 at beginning so energy does not include phi
   local vphi = vpar*Bt/B           --Bt/B = cos(pitch angle). vpar*cos(pitch angle)=vphi
   local toroidal_momentum = m_s*R(r,z)*vphi + q_s*(psi - psi0)

   --Try defining stuff at box minimum
   --local Bm =B_axis*(R_axis/Rmid_min) -- Magnetic field magnitude in the simulation box [T].
   --use actual max
   local Bm = BMax

   local psiShift=0
   local energy = vpar^2/2 + mu*B/m_s
   if (energy - mu*Bm/m_s > 0) then
      psiShift = math.sign(vpar)*m_s*R0*math.sqrt(2*(energy - mu*Bm/m_s))
   end


   --local x = q0/(r0*q_s*B0)*(toroidal_momentum - psiShift)
   local x = 77.8/(q_s)*(toroidal_momentum - psiShift)

   --Calculate the canonical maxwellian
   local n_maxwellian = n0*(0.4*(1.5 + math.tanh(2.*(2.-25.*x)))+0.01)

   local T_maxwellian = T0*((1/3)*(2.+math.tanh(2.*(2.-25.*x)))+0.01)
   return (n_maxwellian/(2*math.pi*T_maxwellian/m_s)^1.5)*math.exp(-m_s*energy/T_maxwellian)
end

function canonical_maxwellian_elc(tcurr,xn)
  return canonical_maxwellian(xn,qe,me,Te0)
end

function canonical_maxwellian_ion(tcurr,xn)
  return canonical_maxwellian(xn,qi,mi,Ti0)
end


-- Source parameters.
-- TRANSP estimated P_OH=1.7919e5 W for this shot.
P_SOL     = 1.7919e5                -- Power into the whole SOL, from experimental heating power [W].
P_src     = P_SOL/ntoroidal         -- Amount of power into flux tube [W].
--n0_src    = 1.                      -- Amplitude of density source [m^{-3} s^{-1}].
n0_src = 4.1e22/5.0
x_src     = xMin                    -- Source start coordinate [m].
sigma_src = 0.2*(x_LCFS-xMin)       -- Characteristic length scale of source [m].
T_srce    = 1.25*Te0                -- Electron source temperature [J].
T_srci    = 1.25*Ti0                -- Ion source temperature [J].

-- Source profiles.
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local srcFloor = 1e-10
   if x < x_src then srcFloor = 1e-2 end -- Higher floor to left of source peak.
   if math.abs(z) < Lz/4 then
      return math.max(n0_src*0.5*(1.-math.tanh(200.*(x-sigma_src))), srcFloor)
   else
      return srcFloor
   end
end
sourceTempIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < x_src + 3*sigma_src then
      return T_srci/5.
   else
      return T_srci*3./8./5.
   end
end
sourceTempElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < x_src + 3*sigma_src then
      --return T_srce
      return T_srce/5.
   else
      --return T_srce*3./8.
      return T_srce*3./8./5.
   end
end

nbc = function(t,xn) return n0 end
Tbc = function(t,xn) return Te0 end
uParbc = function(t,xn) return 0.0 end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

local bcShiftFunc = function(t,xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local r = r_x(x)
   return r0/q0*qprofile(r)*Lz
end

local Nz       = 8
local pOrder   = 1
local dzEff    = Lz/(Nz*(pOrder+1))
local kParMax  = math.pi/dzEff
local kPerpMin = math.min(math.pi/(2.*(xMax-xMin)),math.pi/(2.*Ly))
local omega_H  = kParMax*vte/(kPerpMin*rho_s)
log("\n")
log(string.format("  omega_H = kpar*vte/(kperp*rhos) = %e rad/s\n", omega_H))
log("\n")

local wallTime  = 1.*3600.
local finalTime = .2e-7
local numFrames = 1

plasmaApp = Plasma.App {
   logToFile = true,

   maxWallTime = wallTime,
   tEnd   = finalTime,                -- End time.
   nFrame = numFrames,                -- Number of output frames.
   lower  = {xMin,-Ly/2,-Lz/2},       -- Configuration space lower left.
   upper  = {xMax, Ly/2, Lz/2},       -- Configuration space upper right.
   cells  = {16, 1, Nz},              -- Configuration space cells.
   mapc2p = {grid='grid.bp',mapc2p=function(xc)  -- Transformation from computational to physical coordinates.
      local x, y, z = xc[1], xc[2], xc[3]
      local r = r_x(x)
      -- Map to cylindrical (R, Z, phi) coordinates.
      local R   = R(r, z)
      local Z   = kappa*r*math.sin(z)
      local phi = -q0/r0*y - alpha(r, z, 0)
      -- Map to Cartesian (X, Y, Z) coordinates.
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      return X, Y, Z
   end,},
   basis       = "serendipity",            -- One of "serendipity" or "maximal-order".
   polyOrder   = pOrder,                   -- Polynomial order.
   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.6,
   restartFrameEvery = 0.025/10.0,
   calcIntQuantEvery = 1./(numFrames*10.), -- Aim at 10x more frequently than frames.
   groupDiagnostics  = true, -- False for rescaling

   periodicDirs = {2},     -- Periodic in y only.

   decompCuts = {1,1,1},    -- MPI subdomains/processes.
   parallelizeSpecies = false,

   -- Gyrokinetic electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      lower = {-4.*vte, 0},
      upper = { 4.*vte, me*((4.*vte)^2)/(2.*B0)},
      cells = {8, 4},
      -- Initial conditions.
-- Use local Maxwellian just for the regression test.
--      init = Plasma.FunctionProjection{
----        fromFile = 'elc_restart.in',
--         func = function(t,xn)
--            return canonical_maxwellian(xn, qe,me,Te0)
--         end,
--      },
      init = Plasma.MaxwellianProjection{
         density = function(t, xn) return n0*(0.4*(1.5 + math.tanh(2.*(2.-25.*xn[1])))+0.01) end,
         temperature = function(t,xn) return Te0*((1/3)*(2.+math.tanh(2.*(2.-25.*xn[1])))+0.01) end,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc','ion'},
         frequencies = {nuElc,nuElcIon},
      },
      diff = Plasma.Diffusion {
         coefficient = {function(t,xn) return 1.0 end},
         diffusiveDirs = {1},
         order = 2, -- Diffusion Order: 2,4, or 6.
      },
      --source = Plasma.Source {
      --   kind        = "Maxwellian",  density     = sourceDensity,
      --   --power       = P_src/2,
      --   temperature = sourceTempElc,
      --   diagnostics = {"M0","intM0","M2","intKE"},
      --},
      polarizationDensityFactor = n0,
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "M3par", "M3perp", "intM0", "intM1", "intEnergy"},
      nDistFuncFrame = numFrames/10,
      randomseed = randomseed,

      bcz = {Plasma.AxisymmetricTokamakLimiterBC{xLCFS=x_LCFS,
                diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AxisymmetricTokamakLimiterBC{xLCFS=x_LCFS,
                diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},

-- Use local Maxwellian just for the regression test.
--      bcx = {Plasma.MaxwellianBC{fromFile='elc_bcXlower_Maxwellian_0.in', profile=canonical_maxwellian_elc, kind='canonical',
--      bcx = {Plasma.MaxwellianBC{profile=canonical_maxwellian_elc, kind='canonical',},
--             Plasma.AbsorbBC{}},
      bcx = {Plasma.MaxwellianBC{density=function(t,xn) return n0 end, temperature=function(t,xn) return Te0 end,},
             Plasma.AbsorbBC{}},
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      charge = qi,  mass = mi,
      -- Velocity space grid.
      lower = {-4.*vti, 0},
      upper = { 4.*vti, mi*((4.*vti)^2)/(2.*B0)},
      cells = {8, 4},
      -- Initial conditions.
-- Use local Maxwellian just for the regression test.
--      init = Plasma.FunctionProjection{
----         fromFile='ion_restart.in',
--         func = function(t,xn)
--            return canonical_maxwellian(xn, qi,mi,Ti0)
--         end,
--      },
      init = Plasma.MaxwellianProjection{
         density = function(t, xn) return n0*(0.4*(1.5 + math.tanh(2.*(2.-25.*xn[1])))+0.01) end,
         temperature = function(t,xn) return Ti0*((1/3)*(2.+math.tanh(2.*(2.-25.*xn[1])))+0.01) end,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion','elc'},
         frequencies = {nuIon,nuIonElc},
      },
      diff = Plasma.Diffusion {
         coefficient = {function(t,xn) return 1.0 end},
         diffusiveDirs = {1},
         order = 2, -- Diffusion Order: 2,4, or 6.
      },
      --source = Plasma.Source {
      --   kind        = "Maxwellian",  density     = sourceDensity,
      --   --power       = P_src/2,
      --   temperature = sourceTempIon,
      --   diagnostics = {"M0","intM0","M2","intKE"},
      --},
      polarizationDensityFactor = n0,
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "M3par", "M3perp", "intM0", "intM1", "intEnergy"},
      nDistFuncFrame = numFrames/10,
      randomseed = randomseed,

      bcz = {Plasma.AxisymmetricTokamakLimiterBC{xLCFS=x_LCFS,
                diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AxisymmetricTokamakLimiterBC{xLCFS=x_LCFS,
                diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},

-- Use local Maxwellian just for the regression test.
--      bcx = {Plasma.MaxwellianBC{fromFile='ion_bcXlower_Maxwellian_0.in', profile=canonical_maxwellian_ion, kind='canonical',
--      bcx = {Plasma.MaxwellianBC{profile=canonical_maxwellian_ion, kind='canonical',
--             Plasma.AbsorbBC{}},
      bcx = {Plasma.MaxwellianBC{density=function(t,xn) return n0 end, temperature=function(t,xn) return Ti0 end,},
             Plasma.AbsorbBC{}},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
      -- Dirichlet in x, periodic in y. Potential phi has homogeneous Neumann
      -- BC for the smoothing operation that enforces continuity in z.
      bcLowerPhi = {{T = "N", V = 0.0}, {T = "P"}, {T="AxisymmetricLimitedTokamak", xLCFS=x_LCFS}},
      bcUpperPhi = {{T = "D", V = 0.0}, {T = "P"}, {T="AxisymmetricLimitedTokamak", xLCFS=x_LCFS}},
      isElectromagnetic = false,
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      fromFile='allGeo_restart.in',
      bmag = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         local r = r_x(x)
         local Bt = Bphi(R(r,z))
         local Bp = dPsidr(r,z)/R(r,z)*gradr(r,z)
         return math.sqrt(Bt^2 + Bp^2)
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
