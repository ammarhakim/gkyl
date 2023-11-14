-- Gkyl ------------------------------------------------------------------------
-- 5-moment modeling of dean flow in the 2d axisymmetric (r, z) plane

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Logger = require "Lib.Logger"

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
local T0_mhd = 0.002
-- set one and only one of MA0 and gravity to nonzero
local MA0 = 0  -- vTheta0/vA0; set to nonzero to set a Er x Bz drift
local gravity = 1  -- allows an equilibrium to develop; FIXME unit

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
local Nr, Nz = 60, 300
local decompCuts = {1, 1, 4}
 
local de0 = di0 / math.sqrt(mi/me)
local dz = Lz / Nz
local dr = (r_out - r_inn) / Nr

local tEnd = tA0*0.1
local nFrame = 1

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

momentApp = Moments.App {
   tEnd = tEnd,
   nFrame = nFrame,
   lower = {r_inn, -math.pi/180, 0},
   upper = {r_out, math.pi/180, Lz},
   cells = {Nr, 3, Nz},

   periodicDirs = {3},
   decompCuts = decompCuts,

   mapc2p = function(t, xn)
      local r, th, z = xn[1], xn[2], xn[3]
      return r*math.cos(th), r*math.sin(th), z
   end,

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Moments.Euler { gasGamma = gasGamma },
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
      bcy = { Moments.Field.bcWedge, Moments.Field.bcWedge },
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Moments.Euler { gasGamma = gasGamma },
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
      bcy = { Moments.Field.bcWedge, Moments.Field.bcWedge },
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
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
      is_ext_em_static = true,
      ext_em_func = function (t,xn)
         return 0, 0, 0.0, 0.0, 0.0, Bz0
      end,
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },
      bcy = { Moments.Field.bcWedge, Moments.Field.bcWedge },
   },
}

-- run application
momentApp:run()
