-- Gkyl ------------------------------------------------------------------------
-- 5-moment modeling of dean flow in the 2d axisymmetric (r, z) plane
-- r -> x, theta -> y, z -> z; One cell along theta (y)

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"
local elliptic_integrals = require "elliptic_integrals"

local logger = Logger { logToFile = True }
local log = function(...)
   logger(string.format(...))
   logger("\n")
end

-- basic normalizations
local gasGamma = 1.4
local vA0 = 1.
local mu0 = 1.
local rho0 = 1.
local l0 = 1.
 
local n0 = 1.
local mi = 1.

-- problem-specific mhd parameters
local Bz0_B0 = 1
local a = 1.75  -- Radius of current loops.
local I = 5 * math.pi / mu0 / a  -- Current carried by the loop; Huang2001prl.
local num_extra_coils = 20  -- Number of extra coils above and below the z-ends.
local T0_mhd = 0.002
-- set one and only one of MA0 and gravity to nonzero
local MA0 = 0  -- vTheta0/vA0; set to nonzero to set a Er x Bz drift
local gravity = 0  -- allows an equilibrium to develop; FIXME unit
local gravityDir = 2

local pert = 1e-4
-- math.randomseed(os.time())
math.randomseed(1)

-- non-mhd parameters
local lightSpeed__vA0 = 10
local mi__me = 25.
local Ti0__Te0 = 1.
local di0__l0 = 0.2

-- derived mhd parameters
local B0 = vA0 * math.sqrt(mu0 * rho0)
local Bz0 = B0
local pmag0 = B0*B0 / 2 / mu0
local p0 = n0 * T0_mhd
local tA0 = l0 / vA0
local vTheta0 = vA0 * MA0
 
-- derived non-mhd parameters
local lightSpeed = vA0 * lightSpeed__vA0
local epsilon0 = 1. / lightSpeed^2 / mu0
local me = mi / mi__me
local rhoe0 = rho0 / (1 + mi__me)
local pe0 = p0 / (1 + Ti0__Te0)
local rhoi0 = rho0 - rhoe0
local pi0 = p0 - pe0
local p0 = pe0 + pi0
local rho0 = rhoe0 + rhoi0
local vS0 = math.sqrt(gasGamma*p0/rho0)
local Er0 = -vTheta0 * Bz0

local di0 = l0 * di0__l0
local wpi0 = lightSpeed__vA0 / di0
local qi = wpi0 * math.sqrt(epsilon0 * mi / n0)
local qe = -qi
 
-- domain size
local r_inn = 0.45 * l0
local r_out = 1.45 * l0
local Lz = 5 * l0
local Nr, Nz = 60/2, 300/2
local decompCuts = {1, 1, 2}
 
local de0 = di0 / math.sqrt(mi/me)
local dz = Lz / Nz
local dr = (r_out - r_inn) / Nr

local tEnd = tA0*1
local nFrame = 10

log("%30s = %g", "Lz", Lz)
log("%30s = %g", "r_out-r_inn", r_out-r_inn)
log("%30s = %g", "Lz/de0", Lz/de0)
log("%30s = %g", "(r_out-r_inn)/de0", (r_out-r_inn)/de0)
log("%30s = %g = 1/%g", "dz/de0", dz/de0, de0/dz)
log("%30s = %g = 1/%g", "dr/de0", dr/de0, de0/dr)
log("%30s = %g", "tEnd", tEnd)
log("%30s = %g", "tEnd/tA0", tEnd/tA0)
log("%30s = %g", "nFrame", nFrame)

log("%30s = %g", "vTheta0", vTheta0)
log("%30s = %g", "vTheta0/vA0", vTheta0/vA0)
log("%30s = %g", "vTheta0/vS0", vTheta0/vS0)
log("%30s = %g", "Bz0", Bz0)
log("%30s = %g", "Er0", Er0)

local K = elliptic_integrals.K  -- Complete Elliptic Integral of the First Kind
local E = elliptic_integrals.E  -- Complete Elliptic Integral of the Second Kind

local calcAphiOneCoil = function(r, z, z0)
  --[[Compute A_phi of a current loop using Jackon 3rd ed. eq. (3.37).

  Args:
  r: Number or ndarray. Cylindrical coordinate radius.
  z: Number or ndarray. Cylindrical coordinate z.
  z0: Location of the loop.
  I: Current carried by the loop.
  a: Radius of the loop.

  Returns:
  Aphi: Number of ndarray.
  ]]--
  local zz = z - z0

  local tmp2 = (a + r)^2 + zz^2
  local tmp = math.sqrt(tmp2)
  local k2 = 4*a*r/tmp2
  local k = math.sqrt(k2)

  local Aphi = (mu0*I*a/math.pi) * (1/tmp) * ( (2-k2)*K(k) - 2*E(k) ) / k2
  -- a rearranged form
  -- Aphi = (mu0*I/4/math.pi) * (2/tmp/r) * ( (a^2+r^2+zz^2)*K(k) - tmp^2*E(k) )

  return Aphi
end


local calcAphi = function (r, z)
   local Aphi = 0

   -- Coils; iz=0,1 are the two ends; other iz are extra coils.
   for iz=-num_extra_coils, num_extra_coils+1 do
      Aphi = Aphi + calcAphiOneCoil(r, z, Lz*iz)
   end

   -- Uniform Bz0
   Aphi = Aphi + 0.5*Bz0*r

   return Aphi
end

momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {r_inn, -1, 0},
   upper = {r_out, 1, Lz},
   cells = {Nr, 1, Nz},
   timeStepper = "fvDimSplit",

   periodicDirs = {2, 3},
   decompCuts = decompCuts,

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]
         local rho = rhoe0
         local vr = 0.
         local vt = vTheta0
         local vz = 0.
         local p = pe0 * (1 + pert * math.random())
         local er = p / (gasGamma - 1.) + 0.5 * rho * (vr^2 + vt^2 + vz^2)
         return rho, rho*vr, rho*vt, rho*vz, er
      end,
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]
         local rho = rhoi0
         local vr = 0.
         local vt = vTheta0
         local vz = 0.
         local p = pi0 * (1 + pert * math.random())
         local er = p / (gasGamma - 1.) + 0.5 * rho * (vr^2 + vt^2 + vz^2)
         return rho, rho*vr, rho*vt, rho*vz, er
      end,
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = epsilon0,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]

         local Er = Er0
         local Et = 0
         local Ez = 0
         local Br = 0
         local Bt = 0
         local Bz = 0

         return Er, Et, Ez, Br, Bt, Bz
      end,
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },
   },
 
   axisymmetricFluidSource = Moments.AxisymmetricMomentSource {
      species = {"elc", "ion"},
      timeStepper = "forwardEuler",
      gasGamma = gasGamma,
   },   

   axisymmetricMaxwellSource = Moments.AxisymmetricPhMaxwellSource {
      timeStepper = "forwardEuler",
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "time-centered",
      gravity = gravity,
      gravityDir = gravityDir,
      hasStaticField = true,
      staticEmFunction = function(t, xn)
         local r, z = xn[1], xn[3]
         local Br = -(calcAphi(r, z+0.5*dz) - calcAphi(r, z+0.5*dz)) / dz
         local Bt = 0
         local Bz = (calcAphi(r+0.5*dr, z) - calcAphi(r+0.5*dr, z)) / dr +
                    calcAphi(r, z) / r
         return 0, 0, 0.0, Br, Bt, Bz
      end

   },   
}
-- run application
momentApp:run()
