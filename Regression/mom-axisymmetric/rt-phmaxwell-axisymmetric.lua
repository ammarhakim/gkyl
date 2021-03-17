-- Gkyl ------------------------------------------------------------------------

local Moments = require("App.PlasmaOnCartGrid").Moments()

local mu0 = 1
local epsilon0 = 1
local lightSpeed = 1 / math.sqrt(mu0 * epsilon0)
local chi_e = 0
local chi_m = 1

local l0 = 1
local r_inn = 0.45 * l0
local r_out = 1.45 * l0
local Lz = 5 * l0

local Nr, Nz = 10, 50
local tEnd = 0.1 * l0 / lightSpeed
local nFrame = 1

local momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {r_inn, -1, 0},
   upper = {r_out, 1, Lz},
   cells = {Nr, 1, Nz},
   timeStepper = "fvDimSplit",

   periodicDirs = {2, 3},  -- periodic along theta and z
   decompCuts = {1, 1, 1},

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = epsilon0,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]
         local Er = 0
         local Et = 0
         local Ez = 0
         local Br = 0
         local Bt = 1
         local Bz = 0
         return Er, Et, Ez, Br, Bt, Bz
      end,
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },  -- radial
   },

   axisymSource = Moments.AxisymmetricPhMaxwellSource {
      timeStepper = "forwardEuler",
      epsilon0 = epsilon0,
      mu0 = mu0,
      chi_e = chi_e,
      chi_m = chi_m,
   },   

}

momentApp:run()
