local Moments = require("App.PlasmaOnCartGrid").Moments()

local gasGamma = 3

local Lx = 1
local NX = 2

local cflFrac = 0.9
local tStart = 0.0
local tEnd = 30
local nFrames = 30

local momentApp = Moments.App {

   cflFrac = cflFrac,
   tEnd = tEnd,
   nFrame = nFrames,
   lower = {0},
   upper = {Lx},
   cells = {NX},
   maximumDt = 5e-2,

   periodicDirs = {1},
   decompCuts = nil,

   fluid1 = Moments.Species {
      charge = 1, mass = 5,
      equation = Moments.Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]
         local rho = 2
         local pr = 1 * rho / 5
         local rhou = rho * 1
         local rhov = rho * 0
         local rhow = rho * 0
         local er = pr/(gasGamma-1) + 0.5 * (rhou^2+rhov^2+rhow^2)
         return rho, rhou, rhov, rhow, er
      end,
   },

   fluid2 = Moments.Species {
      charge = 1, mass = 1,
      equation = Moments.Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]
         local rho = 1
         local pr = 2 * rho / 1
         local rhou = rho * 2
         local rhov = rho * 0
         local rhow = rho * 0
         local er = pr/(gasGamma-1) + 0.5 * (rhou^2+rhov^2+rhow^2)
         return rho, rhou, rhov, rhow, er
      end,
   },

}

momentApp:run()
