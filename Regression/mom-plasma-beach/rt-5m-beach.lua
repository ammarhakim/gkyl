-- Gkyl ------------------------------------------------------------------------
local Moments   = require("App.PlasmaOnCartGrid").Moments()
local Euler = Moments.Eq.Euler
local Constants = require "Lib.Constants"

-- Global parameters.
local gasGamma = 5.0/3.0
local xlower= 0.0
local xupper = 1.0
local nx = 400

-- Compute coordinate of interior last edge.
local dx = (xupper-xlower)/nx
local xLastEdge = xupper-dx

local dx100 = (xupper-xlower)/100
-- Compute drive frequency.
local deltaT = dx100/Constants.SPEED_OF_LIGHT
local driveOmega = Constants.PI/10/deltaT

local J0 = 1.0e-12 -- Reference current density (Amps / m^3)

local tEnd = 5.0e-9
local nFrames = 100

epsilon0 = Constants.EPSILON0
mu0 = Constants.MU0

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
   tEnd  = tEnd,
   nFrame = nFrames,
   lower = {xlower},
   upper = {xupper},
   cells = {nx},

   -- Electrons.
   elc = Moments.Species {
      charge = -Constants.ELEMENTARY_CHARGE, mass = Constants.ELECTRON_MASS,

      equation = Euler { gasGamma = gasGamma },
      -- Initial conditions.
      init = function (t, xn)
         local x = xn[1]
         local rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
	    rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
	    Ex, Ey, Ez, Bx, By, Bz = init(x)
         return rho_e, rhovx_e, rhovy_e, rhovz_e, u_e
      end,
      bcx = { Moments.SpeciesBc.bcCopy, Moments.SpeciesBc.bcCopy },
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

      appliedCurrent = function(t, xn)
         local x = xn[1]
	 local J0 = J0
	 local Jy = 0.0
	 if x>xLastEdge then
	    Jy = -J0*math.sin(driveOmega*t)/epsilon0
	 end
	 return 0.0, Jy, 0.0
      end,

      bcx = { Moments.FieldBc.bcCopy, Moments.FieldBc.bcCopy },
   },

}
-- Run application.
momentApp:run()
