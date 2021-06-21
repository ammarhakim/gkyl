-- Gkyl -------------------------------------------------------------------
--
-- Gyrofluid simulation of a single field line with constant
-- density and zero flow velocity but varying parallel temperature.
-- This is a test of the 3+1 Beer+Hammett closure with absorbing BCs.
--
---------------------------------------------------------------------------
local Plasma    = require ("App.PlasmaOnCartGrid").Gyrofluid()
local Constants = require "Lib.Constants"
local Logger    = require "Lib.Logger"

local log = Logger { logToFile = true }

-- Universal constant parameters.
eps0, mu0 = Constants.EPSILON0, Constants.MU0
eV        = Constants.ELEMENTARY_CHARGE
qe, qi    = -eV, eV
me, mp    = Constants.ELECTRON_MASS*16, Constants.PROTON_MASS

-- Plasma parameters.
mi        = 2.014*mp                         -- Deuterium ion mass.
Te0       = 60*eV
n0        = 0.55e17
B0        = 0.242
tau       = 1.
Ti0       = tau*Te0
kperpRhos = 0.3                              -- k_perp*rho_s in the Poisson equation.
beta      = (tau+1)*(2*mu0*n0*Te0)/(B0^2)          -- Plasma beta.

polyOrder = 1
-- Geometry related parameters.
R0 = 1.0
r0 = 0.0
zMin, zMax = -(2.*math.pi*R0)/2., (2.*math.pi*R0)/2.

nuFrac = 1.0  -- Factor multiplying the collisionality.
-- Electron-electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc        = nuFrac*logLambdaElc*(eV^4)*n0/(6*math.sqrt(2)*(math.pi^(3/2))*(eps0^2)*math.sqrt(me)*(Te0^(3/2)))
-- Ion-ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon        = nuFrac*logLambdaIon*(eV^4)*n0/(12*(math.pi^(3/2))*(eps0^2)*math.sqrt(mi)*(Ti0^(3/2)))

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)
c_s = math.sqrt(Te0/mi)

-- Gyrofrequencies and gyroradii.
omega_ci = eV*B0/mi
omega_ce = eV*B0/me
rho_s    = c_s/omega_ci

-- Electron plasma frequency.
omega_pe = math.sqrt((eV^2)*n0/(me*eps0))

-- Perpendicular wavenumber in SI units:
kperp = kperpRhos/rho_s

kpar  = 2.*math.pi/((zMax-zMin))

-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

nu_ee, nu_ii = nuElc, nuIon
kappaParIon  = n0*(vti^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vti*kpar+nu_ii)
kappaPerpIon = n0*(vti^2)/(math.sqrt(2)*D_perp*vti*kpar+nu_ii)
kappaParElc  = n0*(vte^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vte*kpar+nu_ee)
kappaPerpElc = n0*(vte^2)/(math.sqrt(2)*D_perp*vte*kpar+nu_ee)

local Lambda = math.log(mi/(2.*math.pi*me)) -- Sheath factor. Only used by electrons. 

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd       = 0.250e-6,         -- End time.
   nFrame     = 1,               -- Number of output frames.
   lower      = {zMin},          -- Configuration space lower left.
   upper      = {zMax},          -- Configuration space upper right.
   cells      = {64},            -- Configuration space cells.
   basis      = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder  = polyOrder,       -- Polynomial order.
   decompCuts = {1},             -- MPI cuts in each.
   mapc2p     = function(xc)  -- MF (2021/04/20): this is only intended to give a Jacobian=1 for now.
      -- Field-aligned coordinates (x,y).
      local x, y = xc[1], xc[2]
      -- Cylindrical coordinates (R,phi).
      local phi = x/(R0+r0)
      -- Cartesian coordinates (X,Y).
      local X = R0*math.cos(phi)
      local Y = R0*math.sin(phi)
      return X, Y
   end,
   timeStepper = "rk3",              -- One of "rk2" or "rk3".
   cflFrac     = 0.90,
   restartFrameEvery = .05,
   calcIntQuantEvery = 1./10.,  -- Aim at 10x more frequently than frames.

   writeGhost = true,

   -- Gyrofluid ions.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = n0,
         driftSpeed = function(t, zn)
            return math.sqrt((3*Ti0+Te0)/mi)*math.exp(Lambda)*math.tanh(2.*zn[1])
         end,
         perpendicularTemperature = Te0,
         parallelTemperature = function(t, zn)
            return Te0*(1+0.25*math.cos((2.*math.pi/(zMax-zMin))*zn[1]))
         end,
      },
      closure = Plasma.HeatFlux{
         kappaPar = kappaParElc,  kappaPerp = kappaPerpElc,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","upar","Tpar","Tperp","ppar","pperp"},
      bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}}
   },

   -- Gyrofluid electronss.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = n0,
         driftSpeed = function(t, zn)
            return math.sqrt((3*Ti0+Te0)/mi)*math.tanh(2.*zn[1])
         end,
         perpendicularTemperature = Ti0,
         parallelTemperature = function(t, zn)
            return Ti0*(1+0.25*math.cos((2.*math.pi/(zMax-zMin))*zn[1]))
         end,
      },
      closure = Plasma.HeatFlux{
         kappaPar = kappaParIon,  kappaPerp = kappaPerpIon,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","upar","Tpar","Tperp","ppar","pperp"},
      bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}}
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
      kperpSq = kperp^2,
      bcLowerPhi = {{T = "N", V = 0.0}},
      bcUpperPhi = {{T = "N", V = 0.0}},
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function(t,zn) return B0 end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
