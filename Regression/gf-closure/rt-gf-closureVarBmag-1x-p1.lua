-- Gkyl -------------------------------------------------------------------
--
-- Gyrofluid simulation of a single field line with varying B field, constant
-- density and zero flow velocity but varying parallel temperature, and
-- absorbing BCs.
-- This is a test of the 3+1 Beer+Hammett closure.
--
-- The varying magnetic field is like that in a mirror machine.
--
---------------------------------------------------------------------------
local Plasma    = require ("App.PlasmaOnCartGrid").Gyrofluid()
local Constants = require "Lib.Constants"
local Logger    = require "Lib.Logger"
-- SciLua (and extensions) used for AD, root finding and integration.
-- Use sci.math instead of built-in 'math', see https://scilua.org/sci_diff.html.
local math      = require("sci.math").generic
local root      = require "sci.root"
local quad      = require "sci.quad"
local diff      = require "sci.diff-recursive"
local df        = diff.df

local log = Logger { logToFile = true }

polyOrder = 1

-- Universal constant parameters.
eps0 = Constants.EPSILON0
mu0  = Constants.MU0
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV
qi   = eV
me   = Constants.ELECTRON_MASS*16
mp   = Constants.PROTON_MASS

-- Plasma parameters.
mi        = 2.014*mp                         -- Deuterium ion mass.
Te0       = 60*eV
n0        = 0.55e17
B_p       = 0.242
beta      = 0.05                              -- Ratio of plasma to magnetic pressure.
tau       = 1.
Ti0       = tau*Te0
kperpRhos = 0.3                              -- k_perp*rho_s in the Poisson equation.
beta      = (tau+1)*(2*mu0*n0*Te0)/(B_p^2)          -- Plasma beta.

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
omega_ci = eV*B_p/mi
omega_ce = eV*B_p/me
rho_s    = c_s/omega_ci

-- Electron plasma frequency.
omega_pe = math.sqrt((eV^2)*n0/(me*eps0))

-- Geometry parameters.
numCellZ = 64
RatZeq0  = 0.10    -- Radius of the field line at Z=0.
-- Axial coordinate Z extents. Endure that Z=0 is not on
-- the boundary of a cell (due to AD errors).
lowerZ   = -2.501
upperZ   =  2.50
-- Parameters controlling the magnetic equilibrium model.
eqModel  = {
   mcB   = 6.51292/5.,
   gamma = 0.30904,
   Z_m   = 0.98,
}

-- Perpendicular wavenumber in SI units:
kperp = kperpRhos/rho_s

kpar  = 2.*math.pi/(4.*eqModel.Z_m)

-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

nu_ee, nu_ii = nuElc, nuIon
kappaParElc  = n0*(vte^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vte*kpar+nu_ee)
kappaPerpElc = n0*(vte^2)/(math.sqrt(2)*D_perp*vte*kpar+nu_ee)
kappaParIon  = n0*(vti^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vti*kpar+nu_ii)
kappaPerpIon = n0*(vti^2)/(math.sqrt(2)*D_perp*vti*kpar+nu_ii)

local function psi_RZ(RIn,ZIn)
   -- Compute flux function value psi as a function of R and Z.
   -- The input eqModel is a dictionary of model parameters.
   local mcB   = eqModel["mcB"]
   local gamma = eqModel["gamma"]
   local Z_m   = eqModel["Z_m"]

   local psi = 0.5*(RIn^2)*mcB*( 1./(math.pi*gamma*(1.+((ZIn-Z_m)/gamma)^2))
                                +1./(math.pi*gamma*(1.+((ZIn+Z_m)/gamma)^2)) )
   return psi
end

-- Compute the radius R as a function of psi and Z.
-- The input eqModel is a dictionary of model parameters.
local function R_psiZ(psiIn,ZIn)
   local mcB   = eqModel["mcB"]
   local gamma = eqModel["gamma"]
   local Z_m   = eqModel["Z_m"]

   local Rout = math.sqrt(2.*psiIn/(mcB*( 1./(math.pi*gamma*(1.+((ZIn-Z_m)/gamma)^2))
                                         +1./(math.pi*gamma*(1.+((ZIn+Z_m)/gamma)^2)) )))
   return Rout
end

-- Compute the magnetic field (radial, axial and total amplitudes),
-- given a value of the flux function psiIn (which labels the field line)
-- and a 1D array of Z coordinates, ZIn.
local function Bfield_psiZ(psiIn, ZIn)
   local Rcoord = R_psiZ(psiIn,ZIn)

   local mcB   = eqModel["mcB"]
   local gamma = eqModel["gamma"]
   local Z_m   = eqModel["Z_m"]

   local BRad = -(1./2.)*Rcoord*mcB*(-2.*(ZIn-Z_m)/(math.pi*(gamma^3)*((1.+((ZIn-Z_m)/gamma)^2)^2))
                                     -2.*(ZIn+Z_m)/(math.pi*(gamma^3)*((1.+((ZIn+Z_m)/gamma)^2)^2)))

   local BZ   = mcB*( 1./(math.pi*gamma*(1.+((ZIn-Z_m)/gamma)^2))
                     +1./(math.pi*gamma*(1.+((ZIn+Z_m)/gamma)^2)) )

   local Bmag = math.sqrt(BRad^2 + BZ^2)

   return BRad, BZ, Bmag
end

-- Coordinate along the field line (z) as a function of the axial coordinate (Z).
local function z_psiZ(psiIn, ZIn)
   local function integrand(t)
      local _, B_Z, Bmag = Bfield_psiZ(psiIn, t)
      return Bmag/B_Z
   end
   local integral
   local eps = 0.0
   if diff.le(eps, ZIn) then
      integral, _ = quad.dblexp(integrand, eps, ZIn, 1e-14)
   else
      integral, _ = -quad.dblexp(integrand, ZIn, eps, 1e-14)
   end
   return integral
end

-- Stopping criteria for root-finding.
local function stop(tol)
   return function(x, y, xl, xu, yl, yu)
      if diff.lt(math.abs(y), tol) then
         return true
      else
         return false
      end
   end
end

-- Invert z(Z) via root-finding.
local function Z_psiz(psiIn,zIn)
   local function lossF(Z)
      -- Function we will find the roots off. The root corresponds to
      -- the value of Z that satisfies our coordinate transformation
      return zIn - z_psiZ(psiIn,Z)
   end
   local macL = upperZ - lowerZ
   local Zout = 0.0
   local eps  = 1e-14
   if diff.le(eps,zIn) then
      Zout = root.ridders(lossF, eps, upperZ+macL/numCellZ, stop(1e-14))
   else
      Zout = root.ridders(lossF, lowerZ-macL/numCellZ, eps, stop(1e-14))
   end
   return Zout
end

fldLinePsi = psi_RZ(RatZeq0,0.0)         -- Magnetic flux function psi of field line.
zMin       = z_psiZ(fldLinePsi,lowerZ)   -- Min value of coordinate along field line.
zMax       = z_psiZ(fldLinePsi,upperZ)   -- Max value of coordinate along field line.

local nz = 1024
-- Find the throat location with a brute-force search using a finer mesh.
local B_m, Z_m, z_m = -1.e19, 0.0, 0.0
for k = 1, nz+1 do
   local Lz = zMax-zMin              -- Size of the box along the field line.
   local z  = zMin + (k-1)*(Lz/nz)   -- Field-aligned coordinate.
   local Z  = Z_psiz(fldLinePsi,z)   -- Cylindrical axial coordinate Z.
   local R  = R_psiZ(fldLinePsi,Z)   -- Cylindrical radial coordinate R.
   local _, _, Bmag = Bfield_psiZ(fldLinePsi, Z)
   if B_m < Bmag then
      B_m  = Bmag
      Z_m  = math.abs(Z)
      z_m  = math.abs(z)
   end
end

-- Locate the banana tip, assuming particles injected at
-- 45 degree angle at Z=0 (and that all particles have the same
-- pitch angle) so the banatip occurs where B/B(Z=0)=2.
local BmagDiff = 1.e19
local R_bt, Z_bt, z_bt, B_bt = -1.e19, -1.e19, -1.e19, -1.e19
local _, _, BmagAtReq0 = Bfield_psiZ(fldLinePsi, 0.0)
for k = 1, nz+1 do
   local Lp = 2*z_m                  -- Size of the plasma along the field line.
   local z  = -z_m + (k-1)*(Lp/nz)   -- Field-aligned coordinate.
   local Z  = Z_psiz(fldLinePsi,z)   -- Cylindrical axial coordinate Z.
   local R  = R_psiZ(fldLinePsi,Z)   -- Cylindrical radial coordinate R.
   local _, _, Bmag = Bfield_psiZ(fldLinePsi, Z)
   local currBdiff = math.abs(Bmag-2.*BmagAtReq0)
   if currBdiff < BmagDiff then
      BmagDiff = currBdiff
      R_bt = R
      Z_bt = math.abs(Z)
      z_bt = math.abs(z)
      B_bt = Bmag
   end
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd       = 16.0e-8,         -- End time.
   nFrame     = 1,               -- Number of output frames.
   lower      = {zMin},          -- Configuration space lower left.
   upper      = {zMax},          -- Configuration space upper right.
   cells      = {numCellZ},      -- Configuration space cells.
   basis      = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder  = polyOrder,       -- Polynomial order.
   decompCuts = {1},             -- MPI cuts in each.
   -- Dimensions of greater configuration space (for computing metric),
   -- and values at which to evaluate the other coordinates.
   world  = { dim=3, evaluateAt={x=fldLinePsi,y=0.0} },
   mapc2p = function(xc)
      local psi, theta, z = xc[1], xc[2], xc[3]   -- Field-aligned coordinates.

      local Z = Z_psiz(psi,z)   -- Cylindrical axial coordinate Z.
      local R = R_psiZ(psi,Z)   -- Cylindrical radial coordinate R.

      -- Cartesian coordinates on plane perpendicular to Z axis.
      local X = R*math.cos(theta)
      local Y = R*math.sin(theta)

      return X, Y, Z
   end,
   timeStepper = "rk3",              -- One of "rk2" or "rk3".
   cflFrac     = 0.050,
   restartFrameEvery = .05,
   calcIntQuantEvery = 1./10.,  -- Aim at 10x more frequently than frames.

--   periodicDirs = {1},

   -- Gyrofluid ions.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = n0,
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
      bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}},
   },

   -- Gyrofluid electronss.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = n0,
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
      bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}},
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
      bmag = function (t, xn)
         local z   = xn[1]
         local psi = psi_RZ(RatZeq0,0.0)  -- Magnetic flux function psi of field line.

         local Z = Z_psiz(psi, z)
         local _, _, Bmag = Bfield_psiZ(psi, Z)
         return Bmag
      end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
