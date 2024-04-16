-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Mhd = Moments.Eq.Mhd

gasGamma = 5.0/3.0 -- gas adiabatic constant

local function mhd_cons_vars(gasGamma, pv)
   local rho, u, v, w, pr = pv[1], pv[2], pv[3], pv[4], pv[5]
   local q = { }

   q[1] = rho
   q[2] = rho*u
   q[3] = rho*v
   q[4] = rho*w

   q[6] = pv[6]
   q[7] = pv[7]
   q[8] = pv[8]

   local  pb = 0.5*(pv[6]^2 + pv[7]^2 + pv[8]^2)
   q[5] = pr/(gasGamma-1) + 0.5*rho*(u*u+v*v+w*w) + pb

   return q
end
   
-- create app
mhdApp = Moments.App {

   tEnd = 0.5, -- end time
   nFrame = 5, -- number of output frame
   lower = {-0.5, -0.5}, -- lower left corner
   upper = {0.5, 0.5}, -- upper right corner
   cells = {256, 256}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   periodicDirs = { 1, 2 },
   
   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Mhd { 
	 gasGamma = gasGamma,
	 divergenceConstraint = "eight_waves"
      },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local gasGamma = gasGamma
	 local M_PI = math.pi
	 
	 local rho = 25/(36*M_PI)
	 local p = 5/(12*M_PI)
	 local vx = math.sin(2*M_PI*y)
	 local vy = -math.sin(2*M_PI*x)
	 local vz = 0
	 local B0 = 1/math.sqrt(4*M_PI)
	 local Bx = B0*math.sin(2*M_PI*y)
	 local By = B0*math.sin(4*M_PI*x)
	 local Bz = 0
	 local v = {rho, vx, vy, vz, p, Bx, By, Bz }
	 local q = mhd_cons_vars(gasGamma, v)
	 
	 return q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8]
      end,
      evolve = true, -- evolve species?
   },   
}
-- run application
mhdApp:run()
