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
-- Due to axisymmetry, theta direction must have only one cell and be periodic.

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"

local gasGamma = 1.4

local Lx = 1.5
local Ly = 1.0 -- This is arbitrary.
local Lz = 1.0

local NX = 96
local NY = 1  -- This must be 1.
local NZ = 64

local cflFrac = 0.9
local tStart = 0.0
local tEnd = 1.0
local nFrames = 1

local momentApp = Moments.App {
   logToFile = true,

   cflFrac = cflFrac,
   tEnd = tEnd,
   nFrame = nFrames,
   lower = {0, 0, 0},
   upper = {Lx, Ly, Lz},
   cells = {NX, NY, NZ},
   timeStepper = "fvDimSplit",

   periodicDirs = {2},
   decompCuts = decompCuts,

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
      bcx = { Moments.Species.bcAxis, Moments.Species.bcCopy },
      bcz = { Moments.Species.bcWall, Moments.Species.bcWall },
   },

   axisymSource = Moments.AxisymmetricMomentSource {
      species = {"fluid"},
      timeStepper = "semi-exact",
      gasGamma = gasGamma,
   },   
}

momentApp:run()
