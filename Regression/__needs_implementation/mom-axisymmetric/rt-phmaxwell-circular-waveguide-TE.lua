-- 2d axisymmetric (r, z) Maxwell modeling of a circular wave-guide.
-- r -> x, theta -> y, z -> z; One cell along theta (y).
--
-- https://escholarship.org/uc/item/1c49t97t Section 6.1.2
-- https://ntrs.nasa.gov/api/citations/19910007955/downloads/19910007955.pdf Appendix C

local Moments = require("App.PlasmaOnCartGrid").Moments()
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

local logger = Logger { logToFile = True }

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

local pi = math.pi
local xi0 = 3.8317059702075125  -- First zero of J1(x)
local sin = math.sin
local cos = math.cos

local kr = xi0
local kz = pi
local w = math.sqrt(kr^2 + kz^2)
local T = 2*pi/w

local r_inn = 0
local r_out = 1
local Lz = 2*pi/kz
local Nr, Nz = 128, 64
local decompCuts = {1, 1, 1}

local J0 = function (x)
   local ax  = math.abs(x)
   local ans1, ans2, y
   if ax < 8.0 then
      y = x*x
      ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017-y*184.9052456))))
      ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y))))
      return ans1/ans2
   else
      local z = 8.0/ax
      y = z*z
      ans1 = 1.0+y*(-0.1098628627E-2+y*(0.2734510407E-4+y*(-0.2073370639E-5+y*0.2093887211E-6)))
      ans2 = -0.1562499995E-1+y*(0.1430488765E-3+y*(-0.6911147651E-5+y*(0.7621095161E-6-y*0.934935152E-7)))
      local xx = ax-0.785398164
      return math.sqrt(0.636619772/ax)*(math.cos(xx)*ans1-z*math.sin(xx)*ans2)
   end
end

local J1 = function (x)
   local ax = math.abs(x)
   local ans1, ans2, y
   if ax < 8.0 then
      y = x*x
      ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.4826-y*30.16036606)))))
      ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y))))
      return ans1/ans2
   else
      local z = 8.0/ax
      y = z*z
      ans1 = 1.0+y*(0.183105E-2+y*(-0.3516396496E-4+y*(0.2457520174E-5-y*0.240337019E-6)))
      ans2 = 0.04687499995+y*(-0.2002690873E-3+y*(0.8449199096E-5+y*(-0.88228987E-6+y*0.105787412E-6)))
      local xx = ax-2.356194491
      ans1 = math.sqrt(0.636619772/ax)*(math.cos(xx)*ans1-z*math.sin(xx)*ans2)
      return (x >= 0) and ans1 or -ans1
   end
end

local bcRadialEmfAxis = function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
   -- Er=0, Et=0, dEz/dr=0
   qbc[1] = -qin[1]
   qbc[2] = -qin[2]
   qbc[3] = qin[3]

   -- Br=0, Bt=0, d(Bz)/dr=0
   qbc[4] = -qin[4]
   qbc[5] = -qin[5]
   qbc[6] = qin[6]

   qbc[7] = qin[7]
   qbc[8] = qin[8]
end

local bcRadialEmfConductingWall = function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
   local rIn, rOut = xcIn[1], xcOut[1]

   -- Perfectly conducting wall: E_tangent=0, B_normal=0.
   qbc[2] = -qin[2]
   qbc[3] = -qin[3]
   qbc[4] = -qin[4]

   -- Zero surfce charge and current: Er=0, d(r*Bt)/dr=0, d(Bz)/dr=0
   qbc[1] = qin[1]
   qbc[5] = qin[5] * rIn / rOut
   qbc[6] = qin[6]

   qbc[7] = qin[7]
   qbc[8] = qin[8]
end

momentApp = Moments.App {
   logToFile = true,

   tEnd = T,
   nFrame = 4,
   lower = {r_inn, -1, 0},
   upper = {r_out, 1, Lz},
   cells = {Nr, 1, Nz},
   timeStepper = "fvDimSplit",

   periodicDirs = {2, 3},
   decompCuts = decompCuts,

   field = Moments.Field {
      epsilon0 = 1, mu0 = 1,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]
         local tmp = kr/(w^2-kz^2)*J1(kr*r)

         local Br = kz*tmp*sin(kz*z)
         local Bt = 0
         local Bz = J0(kr*r)*cos(kz*z)
         local Er = 0
         local Et = -w*tmp*sin(kz*z)
         local Ez = 0
         return Er, Et, Ez, Br, Bt, Bz
      end,
      bcx = { {bcRadialEmfAxis}, {bcRadialEmfConductingWall} },
   },

   axisymSource = Moments.AxisymmetricPhMaxwellSource {
      timeStepper = "forwardEuler",
   },   
}
momentApp:run()
