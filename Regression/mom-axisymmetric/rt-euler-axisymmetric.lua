-- Axisymmetric modeling of a strong outward-moving shock wave and contact
-- discontinuity and an inward moving rarefaction wave.
--
-- Adapted from https://github.com/ammarhakim/ammar-simjournal/blob/master/sims/s411/s411-riemann-euler-rz.lua
--
-- Originally formulated in Section 3.2 of
-- Langseth, J. O., & LeVeque, R. J. (2000).
-- Journal of Computational Physics, 165(1), 126–166.
--
-- The (r, theta, z) coordinates are treated as (x, y, z).
-- Due to axisymmetry, theta direction must have only one cell and is periodic.

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"

local gasGamma = 1.4

local Lx = 1.5
local Ly = 100.0 --- This is arbitrary.
local Lz = 1.0

local NX = 600
local NY = 1  -- This must be 1.
local NZ = 400
local cfl = 0.9
local tStart = 0.0
local tEnd = 0.7
local nFrames = 100

local momentApp = Moments.App {
   logToFile = true,

   cfl = cfl,
   cflFrac = cflFrac,
   tEnd = tEnd,
   nFrame = nFrames,
   lower = {0, 0, 0},
   upper = {Lx, Ly, Lz},
   cells = {NX, NY, NZ},
   timeStepper = "fvDimSplit",

   periodicDirs = {2},
   decompCuts = {2, 1, 2},

   fluid = Moments.Species {
      charge = 1, mass = 1,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = function (t, xn)
         local r, theta, z = xn[1], xn[2], xn[3]

         local rhoi, pri = 1.0, 5.0
         local rho0, pr0 = 1.0, 1.0
         local rho, pr

         local rloc = math.sqrt(r^2+(z-0.4)^2)
         if (rloc<0.2) then
            rho, pr = rhoi, pri
         else
            rho, pr = rho0, pr0
         end

         return rho, 0.0, 0.0, 0.0, pr/(gasGamma-1)
      end,
      bcx = { 
         {
              BoundaryCondition.Copy { components = {1, 4, 5} },
              BoundaryCondition.Flip { components = {2, 3} },
         },
         Euler.bcCopy
      },
      bcz = { Euler.bcWall, Euler.bcWall },
   },

   axisymSource = Moments.AxisymmetricMomentSource {
      species = {"fluid"},
      timeStepper = "forwardEuler",
      gasGamma = gasGamma,
   },   
}

momentApp:run()
