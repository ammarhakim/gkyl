-- 2d axisymmetric (r, z) perfectly-hyperbolic Maxwell modeling demo.
-- r -> x, theta -> y, z -> z; Theta (y) direction is periodic and has one cell.

local Moments = require("App.PlasmaOnCartGrid").Moments()
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

local logger = Logger { logToFile = True }
local log = function(...)
   logger(string.format(...))
   logger("\n")
end

local mu0 = 1
local epsilon0 = 1

local l0 = 1
local r_inn = 0.45 * l0
local r_out = 1.45 * l0
local Lz = 5 * l0

local lightSpeed = 1 / math.sqrt(epsilon0 * mu0)

local Nr, Nz = 128, 64
local tEnd = 0.1 * l0 / lightSpeed
local nFrame = 1

-- Huang & Hassam (2001). PRL.
local bcRadialEmfConductingWall = function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
   local rIn, rOut = xcIn[1], xcOut[1]

   qbc[1] = qin[1]  -- d(Er)/dr=0
   qbc[2] = -qin[2]  -- Et=0
   qbc[3] = -qin[3]  -- Ez=0

   qbc[4] = -qin[4]  -- Br=0
   qbc[5] = qin[5] * rIn / rOut  -- d(r*Bt)/dr=0
   qbc[6] = qin[5]   -- d(Bz)/dr=0

   -- FIXME Proper BCs for the correction potentials are not well understood.
   qbc[7] = qin[7]
   qbc[8] = qin[8]
end

local momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {r_inn, -1, 0},
   upper = {r_out, 1, Lz},
   cells = {Nr, 1, Nz},
   timeStepper = "fvDimSplit",

   periodicDirs = {2, 3},
   decompCuts = decompCuts,

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]

         local Er = 0
         local Et = 0
         local Ez = 0
         local Br = 0
         local Bt = 0
         local Bz = 0
         return Er, Et, Ez, Br, Bt, Bz
      end,
      bcx = { {bcRadialEmfConductingWall}, {bcRadialEmfConductingWall} },
   },

   axisymSource = Moments.AxisymmetricPhMaxwellSource {
      timeStepper = "forwardEuler",
   },   
}

momentApp:run()
