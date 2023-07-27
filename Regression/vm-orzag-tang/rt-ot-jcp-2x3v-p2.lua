-- Gkyl ------------------------------------------------------------------------
-- 2x3v Orzag Tang simulation from J. Juno, et al. JCP 353, 110 (2018).
--
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = true
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

epsilon0   = 1.0  -- Permittivity of free space.
mu0        = 1.0  -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0)  -- Speed of light.

Te_Ti = 1.0  -- Ratio of electron to ion temperature.
n0    = 1.0  -- Initial number density.

elcMass   = 1.0  -- Electron mass.
elcCharge = -1.0 -- Electron charge.
ionMass   = 25.0 -- Ion mass.
ionCharge = 1.0  -- Ion charge.

vAe   = 0.5
B0    = vAe*math.sqrt(mu0*n0*elcMass)
beta  = 0.08
vtElc = vAe*math.sqrt(beta/2.0)

elcTemp = vtElc^2
ionTemp = elcTemp/Te_Ti

-- Ion velocities
vAi   = vAe/math.sqrt(ionMass)
vtIon = vtElc/math.sqrt(ionMass) --Ti/Te = 1.0

-- Plasma and cyclotron frequencies.
wpe     = math.sqrt(ionCharge^2*n0/(epsilon0*elcMass))
wpi     = math.sqrt(ionCharge^2*n0/(epsilon0*ionMass))
omegaCe = ionCharge*B0/elcMass
omegaCi = ionCharge*B0/ionMass

-- Inertial length.
de = vAe/omegaCe
di = vAi/omegaCi

-- Gyroradii.
rhoi = vtIon/omegaCi
rhoe = vtElc/omegaCe

nuElc = 0.01*omegaCi
nuIon = nuElc/math.sqrt(ionMass)

-- OT initial conditions.
u0x = 0.2*vAi
u0y = 0.2*vAi
B0x = 0.2*B0
B0y = 0.2*B0

-- Domain size and simulation time.
Lx = 20.48*di
Ly = 20.48*di
Nx = 8
Ny = 8
Ncx = 1
Ncy = 1

parSpecies = true
NvElc = 6
NvIon = 6
vMaxElc = 6*vtElc
vMaxIon = 6*vtIon
nFieldFrames = 100
nDistFrames = 100
--tFinal = 100.0/omegaCi
--tFinal = 100*0.0211415 -- 100 steps for 16^2 x 16^3 grid.
tFinal = 40*0.0206722 -- 100 steps for 32^2 x 16^3 grid.
--tFinal = 25*0.0202672 -- 25 steps for 32 x 64 x 16^3 grid.
--tFinal = 25*0.0202672 -- 25 steps for 32 x 64 x 16^3 grid.

-- Estimate dt and number of steps.
dx = Lx/Nx
polyOrder = 2
deltaT = 2*(dx/lightSpeed)/(2*polyOrder+1)
nSteps = tFinal/deltaT

log("%50s = %g", "mi/me", ionMass / elcMass)
log("%50s = %g", "wpe/OmegaCe", wpe / omegaCe)
log("%50s = %g", "proton beta", beta)
log("%50s = %g", "vte/c", vtElc / lightSpeed)
log("%50s = %g", "vti/c", vtIon / lightSpeed)
log("%50s = %g", "electron plasma frequency (wpe) ", wpe)
log("%50s = %g", "electron cyclotron frequency (OmegaCe) ", omegaCe)
log("%50s = %g", "ion plasma frequency (wpi) ", wpi)
log("%50s = %g", "ion cyclotron frequency (OmegaCi) ", omegaCi)
log("%50s = %g", "electron-electron collision frequency (nuee) ", nuElc)
log("%50s = %g", "ion-ion collision frequency (nuii) ", nuIon)
log("%50s = %g", "electron inertial length (de) ", de)
log("%50s = %g", "ion inertial length (di) ", di)
log("%50s = %g", "Max ion velocity/vti", vMaxIon/vtIon)
log("%50s = %g", "Max elc velocity/vte", vMaxElc/vtElc)
log("%50s = %g", "Number of grid cells per di in x", Nx/(Lx/di))
log("%50s = %g", "Number of grid cells per de in x", Nx/(Lx/de))
log("%50s = %g", "tFinal ", tFinal)
log("%50s = %g", "End time in inverse ion cyclotron periods", tFinal*omegaCi)
log("%50s = %g", "Estimated time step", deltaT)
log("%50s = %g", "Estimated number of time steps", nSteps)

-- Maxwellian in 2x3v.
local function maxwellian3D(n, vx, vy, vz, ux, uy, uz, temp, mass)
   local v2 = (vx - ux)^2 + (vy - uy)^2 + (vz - uz)^2
   return n/math.sqrt((2*math.pi*temp/mass)^3)*math.exp(-mass*v2/(2.0*temp))
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd        = tFinal,        -- End time.
   nFrame      = nFieldFrames,  -- Number of output frames.
   lower       = {0.0, 0.0},    -- Configuration space lower left.
   upper       = {Lx, Ly},      -- Configuration space upper right.
   cells       = {Nx, Ny},      -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,     -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2" or "rk3".

   -- Parallelization:
   decompCuts         = {Ncx, Ncy},  -- MPI domains in each configuration direction.
   parallelizeSpecies = parSpecies,

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2},  -- Periodic directions.
   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.001,
   -- Frequency with which to write restart frames.
   restartFrameEvery = 0.01,

   -- Electrons.
   elc = Plasma.Species {
      charge = elcCharge,  mass = elcMass,
      -- Velocity space grid.
      lower = {-vMaxElc, -vMaxElc, -vMaxElc},
      upper = { vMaxElc,  vMaxElc,  vMaxElc},
      cells = {NvElc, NvElc, NvElc},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0
   
         local vdrift_x = -u0x*sin(_2pi*y/Ly)
         local vdrift_y = u0y*sin(_2pi*x/Lx)
         local vdrift_z = -Jz / qi

         local fv = maxwellian3D(n0, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, elcTemp, elcMass)
         
         return fv
      end,
      nDistFuncFrame = nDistFrames,
      evolve = true,  -- Evolve species?
--      diagnostics = { "M0", "M1i", "M2", "M2ij", "M3i", "intM0", "intM1i", "intM2Flow", "intM2Thermal" }, 
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },

   -- Protons.
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- Velocity space grid
      lower = {-vMaxIon, -vMaxIon, -vMaxIon},
      upper = { vMaxIon,  vMaxIon,  vMaxIon},
      cells = {NvIon, NvIon, NvIon},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local vdrift_x = -u0x*sin(_2pi*y/Ly)
         local vdrift_y = u0y*sin(_2pi*x/Lx)
         local vdrift_z = 0.0

         local fv = maxwellian3D(n0, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, ionTemp, ionMass)

         return fv
      end,
      nDistFuncFrame = nDistFrames,
      evolve = true,  -- Evolve species?
--      diagnostics = { "M0", "M1i", "M2", "M2ij", "M3i", "intM0", "intM1i", "intM2Flow", "intM2Thermal" }, 
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0

         local Bx = -B0x*sin(_2pi*y/Ly)
         local By = B0y*sin(_4pi*x/Lx)
         local Bz = B0

         -- Assumes qi = abs(qe)
         local u_xe = -u0x*sin(_2pi*y/Ly)
         local u_ye = u0y*sin(_2pi*x/Lx)
         local u_ze = -Jz/qi
   
         -- E = - v_e x B ~  (J - u) x B
         local Ex = - (u_ye*Bz - u_ze*By)
         local Ey = - (u_ze*Bx - u_xe*Bz)
         local Ez = - (u_xe*By - u_ye*Bx)

         return Ex, Ey, Ez, Bx, By, Bz
      end,

      evolve = true,  -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
