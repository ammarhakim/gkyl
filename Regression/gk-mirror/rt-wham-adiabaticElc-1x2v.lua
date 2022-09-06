-- Gkyl ------------------------------------------------------------------------
--
-- 1x2v gyrokinetic calculation in a mirror geometry modeled after WHAM.
-- Notes: 
--   a) Because of some limitations of AD in the calculation of the
--      Jacobian near Z=0 we do not center the domain at Z=0.
--   b) In order to compute the right metric quantities, we must do so
--      in 3D and then evaluate them at a value in the other two coordinates.
--
-- All dimensional quantities are in SI units unless noted otherwise.
--
-- Manaure Francisquez.
-- January 2020.
--
--------------------------------------------------------------------------------

local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic()
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

log("\n")
log("  1x2v gyrokinetic WHAM simulation\n")
log("\n")

polyOrder = 1

-- Universal constant parameters.
eps0, mu0 = Constants.EPSILON0, Constants.MU0
eV        = Constants.ELEMENTARY_CHARGE
qe, qi    = -eV, eV
me, mp    = Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters.
mi        = 2.014*mp                         -- Deuterium ion mass.
Te0       = 940*eV
n0        = 3.e19
B_p       = 0.53
beta      = 0.4                              -- Ratio of plasma to magnetic pressure.
tau       = (B_p^2)*beta/(2*mu0*n0*Te0)-1    -- Ti/Te ratio.
Ti0       = tau*Te0
kperpRhos = 0.3                              -- k_perp*rho_s in the Poisson equation.

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

log(string.format("  Collision frequencies (including nuFrac=%f):\n",nuFrac))
log(string.format("    nu_ee = %e Hz\n", nuElc))
log(string.format("    nu_ii = %e Hz\n", nuIon))
log("\n")

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)
c_s = math.sqrt(Te0/mi)

-- Gyrofrequencies and gyroradii.
omega_ci = eV*B_p/mi
rho_s    = c_s/omega_ci 

-- Perpendicular wavenumber in SI units:
kperp = kperpRhos/rho_s

-- Geometry parameters.
numCellZ = 128
RatZeq0  = 0.10    -- Radius of the field line at Z=0.
-- Axial coordinate Z extents. Endure that Z=0 is not on
-- the boundary of a cell (due to AD errors).
lowerZ   = -2.501
upperZ   =  2.50
-- Parameters controlling the magnetic equilibrium model.
eqModel  = {
   mcB   = 6.51292,
   gamma = 0.124904,
   Z_m   = 0.98,
}

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

-- Print banana tip location and magnetic field there.
log("  Banana tip field & location:\n")
log(string.format("    B_bt = %f T\n", B_bt))
log(string.format("    R_bt = %f m\n", R_bt))
log(string.format("    Z_bt = %f m\n", Z_bt))
log(string.format("    z_bt = %f m\n", z_bt))


-- Some physics parameters at mirror throats.
local n_m, Te_m = 0.0, 0.0
local nz = 1024
for k = 1, nz+1 do
   local Lz = zMax-zMin              -- Size of the box along the field line.
   local z  = zMin + (k-1)*(Lz/nz)   -- Field-aligned coordinate.
   local Z  = Z_psiz(fldLinePsi,z)   -- Cylindrical axial coordinate Z.
   local R  = R_psiZ(fldLinePsi,Z)   -- Cylindrical radial coordinate R.
   local _, _, Bmag = Bfield_psiZ(fldLinePsi, Z)
   if math.abs(B_m-Bmag) < 1.e-10 then
      n_m  = n0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[2])
      Te_m = Te0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[2])
   end
end
local Ti_m = tau*Te_m
local cs_m = math.sqrt(Te_m*(1.+tau)/mi)

-- Print mirror throat location and magnetic field there.
log("  At mirror throat:\n")
log(string.format("    B_m  = %f T\n", B_m))
log(string.format("    Z_m  = %f m\n", Z_m))
log(string.format("    z_m  = %f m\n", z_m))
log(string.format("    n_m  = %e m^{-3}\n", n_m))
log(string.format("    Te_m = %f eV\n", Te_m/eV))
log(string.format("    Ti_m = %f eV\n", Ti_m/eV))
log(string.format("    cs_m = %e m/s\n", cs_m))

-- Source parameters.
NSrcIon      = 3.1715e23/8.
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
   local z   = xn[1]
   local psi = psi_RZ(RatZeq0,0.0)   -- Magnetic flux function psi of field line.
   local Z   = Z_psiz(psi,z)         -- Cylindrical axial coordinate.

   local NSrc, zSrc        = NSrcElc, zSrcElc   
   local sigSrc, NSrcFloor = sigSrcElc, NSrcFloorElc
   if math.abs(Z) <= Z_m then
      return math.max(NSrcFloor,(NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2))))
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
   local z   = xn[1]
   local psi = psi_RZ(RatZeq0,0.0)   -- Magnetic flux function psi of field line.
   local Z   = Z_psiz(psi,z)         -- Cylindrical axial coordinate.

   local NSrc, zSrc        = NSrcIon, zSrcIon   
   local sigSrc, NSrcFloor = sigSrcIon, NSrcFloorIon
   if math.abs(Z) <= Z_m then
      return math.max(NSrcFloor,(NSrc/math.sqrt(2.*math.pi*(sigSrc^2)))*math.exp(-((z-zSrc)^2)/(2.*(sigSrc^2))))
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

local dzEff   = (zMax-zMin)/(numCellZ*(polyOrder+1))
local kParMax = math.pi/dzEff
local omega_H = kParMax*vte/kperpRhos

log("\n")
log(string.format("  omega_H = kpar*vte/(kperp*rhos) = %e rad/s\n", omega_H))
log("\n")

-- Simulation App.
plasmaApp = Plasma.App {
   logToFile = true,

   tEnd       = 0.50e-8,                   -- End time.
   nFrame     = 1,                      -- Number of output frames.
   lower      = {zMin},                  -- Configuration space lower left.
   upper      = {zMax},                  -- Configuration space upper right.
   cells      = {numCellZ},              -- Configuration space cells.
   basis      = "serendipity",           -- One of "serendipity" or "maximal-order".
   polyOrder  = polyOrder,               -- Polynomial order.
   decompCuts = {1},
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

   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.8,
--   restartFrameEvery = .05,
--   calcIntQuantEvery = 1./400.,

   -- Adiabatic electrons.
   elc = Plasma.AdiabaticSpecies {
      charge = qe, mass = me,
      temp   = Te0,
      init   = function (t, xn) return n0 end,  -- Initial conditions. Not used with ASheathPotential.
      evolve = false, -- Evolve species?
   },

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      lower = {-3.75*vti, 0},
      upper = { 3.75*vti, mi*((3.*vti)^2)/(2*B_p)},
      cells = {64, 8},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
--         fromFile = "ionIC.bp",
         density = function (t, xn)
            local z   = xn[1]
            local psi = psi_RZ(RatZeq0,0.0)   -- Magnetic flux function psi of field line.
            local Z   = Z_psiz(psi,z)         -- Cylindrical axial coordinate.
            local R   = R_psiZ(psi,Z)         -- Cylindrical radial coordinate.
            local _, _, Bmag = Bfield_psiZ(psi, Z)
            if math.abs(Z) <= Z_bt then
               return n0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return n0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[2])
            else
               return n_m*math.sqrt(Bmag/B_m)
            end
         end,
         driftSpeed = function (t, xn)
            local z = xn[1]
            return cs_m*(z/z_m)
         end,
         temperature = function (t, xn)
            local z   = xn[1]
            local psi = psi_RZ(RatZeq0,0.0)   -- Magnetic flux function psi of field line.
            local Z   = Z_psiz(psi,z)         -- Cylindrical axial coordinate.
            local R   = R_psiZ(psi,Z)         -- Cylindrical radial coordinate.
            local _, _, Bmag = Bfield_psiZ(psi, Z)
            if math.abs(Z) <= Z_bt then
               return Ti0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return Ti0*math.sqrt((1. - ((R-R_bt)/alim)^2)^alphaIC[2])
            else
               return Ti_m*math.sqrt(Bmag/B_m)
            end
         end,
         scaleWithSourcePower = false,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
      source = Plasma.Source {
--         fromFile    = "ion_fSourceIC.bp",
         density     = srcDenIon,
         temperature = srcTempIon,
         diagnostics = {"twoFiles","intM0", "intM2"},
      },
      evolve = true, -- Evolve species?
      diagnostics = {"twoFiles","M0", "Upar", "Temp", "Tperp", "Tpar", "intM0", "intM1", "intM2", "intKE", "intEnergy" },
      randomseed = randomseed,
      bcx = {Plasma.SheathBC{diagnostics={"twoFiles","M0","M1","M2","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.SheathBC{diagnostics={"twoFiles","M0","M1","M2","Energy","intM0","intM1","intKE","intEnergy"}}},
   },

   -- Field solver.
   field = Plasma.AmbipolarSheathField {
      bcLowerPhi = {{ T ="N", V = 0.0}},
      bcUpperPhi = {{ T ="N", V = 0.0}},
      evolve     = true,   -- Evolve fields?
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local z   = xn[1]
         local psi = psi_RZ(RatZeq0,0.0)  -- Magnetic flux function psi of field line.

         local Z = Z_psiz(psi, z)
         local _, _, Bmag = Bfield_psiZ(psi, Z)
         return Bmag
      end,
      evolve   = false,   -- Geometry is not time-dependent.
--      fromFile = "allGeoIC.bp",
   },

}
-- Run application.
plasmaApp:run()
