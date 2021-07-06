-- Gkyl -------------------------------------------------------------------
--
-- Gyrofluid simulation of a single open field line with sheath BCs.
--
---------------------------------------------------------------------------
local Plasma    = require ("App.PlasmaOnCartGrid").Gyrofluid()
local Constants = require "Lib.Constants"
local Logger    = require "Lib.Logger"

polyOrder = 1

-- Universal constant parameters.
eps0, mu0 = Constants.EPSILON0, Constants.MU0
eV = Constants.ELEMENTARY_CHARGE
qe, qi = -eV, eV
me, mp = 16*Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters.
mi        = 2.014*mp                         -- Deuterium ion mass.
Te0       = 60*eV
n0        = 0.55e17
B0        = 0.242
beta      = 0.05                              -- Ratio of plasma to magnetic pressure.
tau       = 1.
Ti0       = tau*Te0
kperpRhos = 0.3                              -- k_perp*rho_s in the Poisson equation.
beta      = (tau+1)*(2*mu0*n0*Te0)/(B0^2)          -- Ti/Te ratio.

-- Parameters controlling initial conditions.
alim    = 0.125    -- Plasma limiting radius.
alphaIC = {2,10}

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

-- Geometry parameters.
numCellZ = 64
R0 = 1.0
r0 = 0.0
zMin = -2.50   -- Min value of coordinate along field line.
zMax =  2.50   -- Max value of coordinate along field line.
Z_m  = zMax/2.

-- Perpendicular wavenumber in SI units:
kperp = kperpRhos/rho_s

kpar  = 2.*math.pi/(4.*Z_m)
-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

nu_ee, nu_ii = nuElc, nuIon
kappaParIon  = n0*(vti^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vti*kpar+nu_ii)
kappaPerpIon = n0*(vti^2)/(math.sqrt(2)*D_perp*vti*kpar+nu_ii)
kappaParElc  = n0*(vte^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vte*kpar+nu_ee)
kappaPerpElc = n0*(vte^2)/(math.sqrt(2)*D_perp*vte*kpar+nu_ee)

-- Source parameters.
NSrcIon      = 3.1715e23/80.
zSrcIon      = 0.0
sigSrcIon    = Z_m/4.
NSrcFloorIon = 0.05*NSrcIon
TSrc0Ion     = Ti0*1.25
TSrcFloorIon = TSrc0Ion/8.
NSrcElc      = NSrcIon
zSrcElc      = zSrcIon
sigSrcElc    = sigSrcIon
NSrcFloorElc = NSrcFloorIon
TSrc0Elc     = TSrc0Ion/tau
TSrcFloorElc = TSrcFloorIon/tau

-- Source functions.
srcDenElc = function(t, xn)
   local z = xn[1]
   local NSrc, zSrc        = NSrcElc, zSrcElc
   local sigSrc, NSrcFloor = sigSrcElc, NSrcFloorElc
   if math.abs(z) <= Z_m then
--      return math.max(NSrcFloor,(NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2))))
      return (NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2)))
   else
      return 1.e-16
   end
end
srcTempElc = function(t, xn)
   local z = xn[1]
   local sigSrc = sigSrcElc
   local TSrc0, Tfloor = TSrc0Elc, TSrcFloorElc
   if math.abs(z) <= 2.*sigSrc then
      return TSrc0
   else
      return Tfloor
   end
end
srcDenIon = function(t, xn)
   local z = xn[1]
   local NSrc, zSrc        = NSrcIon, zSrcIon
   local sigSrc, NSrcFloor = sigSrcIon, NSrcFloorIon
   if math.abs(z) <= Z_m then
--      return math.max(NSrcFloor,(NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2))))
      return (NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2)))
   else
      return 1.e-16
   end
end
srcTempIon = function(t, xn)
   local z = xn[1]
   local sigSrc = sigSrcIon
   local TSrc0, Tfloor = TSrc0Ion, TSrcFloorIon
   if math.abs(z) <= 2.*sigSrc then
      return TSrc0
   else
      return Tfloor
   end
end

local Lambda = math.log(mi/(2.*math.pi*me)) -- Sheath factor. Only used by electrons.

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd       = 6.0e-6,         -- End time.
   nFrame     = 1,              -- Number of output frames.
   lower      = {zMin},          -- Configuration space lower left.
   upper      = {zMax},          -- Configuration space upper right.
   cells      = {numCellZ},      -- Configuration space cells.
   basis      = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder  = polyOrder,       -- Polynomial order.
   decompCuts = {1},             -- MPI cuts in each.
   mapc2p = function(xc)
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
   cflFrac     = 0.50,
   restartFrameEvery = .05,
   calcIntQuantEvery = 1./10.,  -- Aim at 10x more frequently than frames.
--   writeGhost = true,

   -- Gyrofluid ions.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = function (t, xn)
            local z = xn[1]
            return math.max(n0*1.e-8,(n0/math.sqrt(2.*math.pi*(sigSrcIon^2)))*math.exp(-((z-0.)^2)/(2.*(sigSrcIon^2))))
         end,
         driftSpeed = function (t, xn)
            local z = xn[1]
            return math.sqrt((3.*Ti0+Te0)/mi)*math.exp(Lambda)*(z/zMax)
         end,
         perpendicularTemperature = Te0,
         parallelTemperature = Te0,
      },
--      closure = Plasma.HeatFlux{
--         kappaPar = kappaParElc,  kappaPerp = kappaPerpElc,
--      },
--      coll = Plasma.PASCollisions {
--         collideWith = {'elc'},
--         frequencies = {nuElc},
--      },
      source = Plasma.Source {
--         fromFile    = "ion_fSourceIC.bp",
         density                  = srcDenElc,
         parallelTemperature      = srcTempElc,
         perpendicularTemperature = srcTempElc,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","Upar","Tpar","Tperp","Ppar","Pperp"},
      bcx = {Plasma.SheathBC{}, Plasma.SheathBC{}},
   },

   -- Gyrofluid electronss.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = function (t, xn)
            local z = xn[1]
            return math.max(n0*1.e-8,(n0/math.sqrt(2.*math.pi*(sigSrcElc^2)))*math.exp(-((z-0.)^2)/(2.*(sigSrcElc^2))))
         end,
         driftSpeed = function (t, xn)
            local z = xn[1]
            return math.sqrt((3.*Ti0+Te0)/mi)*(z/zMax)
         end,
         perpendicularTemperature = Ti0,
         parallelTemperature = Ti0,
      },
      source = Plasma.Source {
--         fromFile    = "ion_fSourceIC.bp",
         density                  = srcDenIon,
         parallelTemperature      = srcTempIon,
         perpendicularTemperature = srcTempIon,
      },
--      closure = Plasma.HeatFlux{
--         kappaPar = kappaParIon,  kappaPerp = kappaPerpIon,
--      },
--      coll = Plasma.PASCollisions {
--         collideWith = {'ion'},
--         frequencies = {nuIon},
--      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","Upar","Tpar","Tperp","Ppar","Pperp"},
      bcx = {Plasma.SheathBC{}, Plasma.SheathBC{}},
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
--      fromFile = "allGeo_restart.bp",
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
