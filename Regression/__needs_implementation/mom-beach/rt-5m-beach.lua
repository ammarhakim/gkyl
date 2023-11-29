-- Gkyl ------------------------------------------------------------------------
local Moments   = require("App.PlasmaOnCartGrid").Moments()
local Constants = require "Lib.Constants"

-- Global parameters.
local cfl      = 0.4
local gasGamma = 5.0/3.0
local xlower   = 0.0
local xupper   = 1.0
local nx       = 400

-- Compute coordinate of interior last edge.
local dx        = (xupper-xlower)/nx
local xLastEdge = xupper-dx

local dx100 = (xupper-xlower)/100
-- Compute drive frequency.
local deltaT     = dx100/Constants.SPEED_OF_LIGHT
local driveOmega = Constants.PI/10/deltaT

local tEnd    = 5.0e-9
local nFrames = 100

epsilon0 = Constants.EPSILON0
mu0      = Constants.MU0

local function init(x)
   local wpdt = 25*(1-x)^5   -- Plasma frequency.
   local factor = deltaT^2*Constants.ELEMENTARY_CHARGE^2/(Constants.ELECTRON_MASS*Constants.EPSILON0)
   local ne = wpdt^2/factor
   local te = 1.0*Constants.EV2KELVIN   -- Electron temperature [K].
   local pre = ne*Constants.BOLTZMANN_CONSTANT*te
   return Constants.ELECTRON_MASS*ne, 0, 0, 0, pre/(gasGamma-1),
          Constants.PROTON_MASS*ne, 0, 0, 0, pre/(gasGamma-1), 
          0, 0, 0, 0, 0, 0, 0, 0
end

momentApp = Moments.App {
   tEnd   = tEnd,
   nFrame = nFrames,
   lower  = {xlower},
   upper  = {xupper},
   cells  = {nx},
   cfl    = cfl,
   maximumDt   = 3e-12,

   -- Decomposition for configuration space.
   decompCuts = {1},    -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {},   -- Periodic directions.

   -- Electrons.
   elc = Moments.Species {
      charge = -Constants.ELEMENTARY_CHARGE, mass = Constants.ELECTRON_MASS,

      equation    = Moments.Euler { gasGamma = gasGamma },
      -- Initial conditions.
      init = function (t, xn)
         local x = xn[1]
         local rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
               rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
               Ex, Ey, Ez, Bx, By, Bz = init(x)
         return rho_e, rhovx_e, rhovy_e, rhovz_e, u_e
      end,
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },


   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x = xn[1]
         local rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
               rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
               Ex, Ey, Ez, Bx, By, Bz = init(x)
         return Ex, Ey, Ez, Bx, By, Bz
      end,
      bcx    = { Moments.Field.bcCopy, Moments.Field.bcCopy },
   },

   -- emSource = Moments.CollisionlessEmSource {
   --    species = {"elc"},
   --    timeStepper = "time-centered",
   --
   --    -- Additional source terms.
   --    -- Not enabled here; for demo purpose.
   --    -- Note: Dp are c arrays and use 0-based indices; xc and qbym are lua
   --    --       arrays and use 1-based indices
   --    hasAuxSourceFunction = true,
   --    auxSourceFunction    = function (
   --       self, xc, t, epsilon0, qbym, fDp, emDp, auxSrcDp)
   --       local x = xc[1]
   --       local nFluids = #qbym
   --
   --       local J0 = 1.0e-12   -- Amps/m^3.
   --       -- Auxiliary source for currents.
   --       for s=0, 0 do
   --          auxSrcDp[s*3+0] = 0
   --          auxSrcDp[s*3+1] = 0
   --          auxSrcDp[s*3+2] = 0
   --       end
   --
   --       -- Auxiliary source for E field.
   --       auxSrcDp[nFluids*3+0] = 0
   --       auxSrcDp[nFluids*3+1] = 0
   --       auxSrcDp[nFluids*3+2] = 0
   --
   --       if (x>xLastEdge) then
   --          auxSrcDp[nFluids*3+1] = -J0*math.sin(driveOmega*t)/Constants.EPSILON0
   --       end
   --    end,
   -- },   

}
-- Run application.
momentApp:run()
