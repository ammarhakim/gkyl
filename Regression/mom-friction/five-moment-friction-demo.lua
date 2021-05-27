local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"

local gasGamma = 3

local Lx = 1
local NX = 2

local cflFrac = 0.9
local tStart = 0.0
local tEnd = 30
local nFrames = 30

local momentApp = Moments.App {
   logToFile = true,

   cflFrac = cflFrac,
   tEnd = tEnd,
   nFrame = nFrames,
   lower = {0},
   upper = {Lx},
   cells = {NX},
   timeStepper = "fvDimSplit",
   maximumDt = 5e-2,

   periodicDirs = {1},
   decompCuts = nil,

   fluid1 = Moments.Species {
      charge = 1, mass = 5,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
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
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
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

   friction = Moments.MomentFrictionSource {
      species = {"fluid1", "fluid2"},
      timeStepper = "time-centered", -- time-centered, exact, or forwardEuler
      gasGamma = gasGamma,
      nu = {0.1},  -- In the order like nu12, nu13, nu14, nu23, nu24, nu34.
                   -- Other values like nu21 will be computed during the
                   -- simulation to satisfy rho1*nu12=rho2*nu21 due to momentum
                   -- conservation.
      updatePressure = true,
   },   
}

momentApp:run()
