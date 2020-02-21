-- Gkyl ------------------------------------------------------------------------

local Proto = require("Proto")
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local Constants = require "Constants"
local Logger = require "Logger"

local logger = Logger {
   logToFile = True
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

----------------
-- PARAMETERS --
----------------
-- physical constants
               gasGamma = 5/3
             lightSpeed = Constants.SPEED_OF_LIGHT/100
                    mu0 = Constants.MU0
               epsilon0 = 1/mu0/(lightSpeed^2)
-- planet radius
                     Re = 2634.1e3
-- planet dipole streghth
                     Dx =    -18e-9*Re^3
                     Dy =   51.8e-9*Re^3
                     Dz = -716.8e-9*Re^3
                     r1 = 0.5*Re -- zero B field when r < r1
-- radius of inner boundary
                     r0 = 1*Re
-- inflows parameters
                 rho_in = 56e6 * Constants.PROTON_MASS
                  vx_in = 140e3
                  vy_in = 0
                  vz_in = 0
                   p_in = 3.8e-9
                  Bx_in = 0
                  By_in = -6e-9
                  Bz_in = -77e-9
-- mirdip setup specfications
                r_ramp1 = 2*Re
                r_ramp2 = 2.5*Re
                stretch = 1
                   xmir = -2.578125*Re
-- kinetic parameters
              numFluids = 2
             mass_ratio = 25
         pressure_ratio = 5
                 d_i_in = 0.2*Re
                ionMass = 14*Constants.PROTON_MASS
                elcMass = ionMass/mass_ratio
              ionCharge = ionMass / d_i_in / math.sqrt(mu0*rho_in)
              elcCharge = -ionCharge
                 charge = {elcCharge, ionCharge}
                   mass = {elcMass, ionMass}
                 d_e_in = d_i_in * math.sqrt(elcMass / ionMass)
-- domain and grid
          numDimensions = 3
           xlo, xup, nx = -24*Re, 24*Re, 96
           ylo, yup, ny = -24*Re, 24*Re, 96
           zlo, zup, nz = -24*Re, 24*Re, 96
                  lower = {xlo, ylo, zlo}
                  upper = {xup, yup, zup}
                  cells = {nx, ny, nz}
-- nonuniform grid
-- f_dx is inverse of a two-peak distribution.
-- x1, wx1, and rx1 are the center, Gaussian width, and max/min ratio
-- of the first Gaussian-like distribution.
-- f_dy and f_dz are single-peak distributions.
-- wy, fwy, ry are the Gaussian-width, flat region width, and max/min
-- ratio of the distribution.
      useNonUniformGrid = false
              wx,fwx,rx =   1.9*Re, 2.25*Re, 50
              wy,fwy,ry =    2.*Re,   2.*Re, 50
              wz,fwz,rz =    2.*Re,   2.*Re, 50

-- computational constants
                    cfl = 0.9
                   cflm = 1.1*cfl
                limiter = "monotonized-centered"
-- i/o control
                   tEnd = 2700
                 tFrame = 60
                 nFrame = math.floor(tEnd/tFrame)
-- derived parameters
               rho_e_in = rho_in / (1 + mass_ratio)
             rhovx_e_in = rho_e_in * vx_in
             rhovy_e_in = rho_e_in * vy_in
             rhovz_e_in = rho_e_in * vz_in
                 p_e_in = p_in / (1 + pressure_ratio)
                 e_e_in = p_e_in / (gasGamma - 1.) + 0.5 * (rhovx_e_in^2 + rhovy_e_in^2 + rhovz_e_in^2) / rho_e_in
               Pxx_e_in = p_e_in + vx_in^2 * rho_e_in
               Pyy_e_in = p_e_in + vy_in^2 * rho_e_in
               Pzz_e_in = p_e_in + vz_in^2 * rho_e_in
               Pxy_e_in = vx_in * vy_in * rho_e_in
               Pxz_e_in = vx_in * vz_in * rho_e_in
               Pyz_e_in = vy_in * vz_in * rho_e_in

               rho_i_in = rho_in - rho_e_in
             rhovx_i_in = rho_i_in * vx_in
             rhovy_i_in = rho_i_in * vy_in
             rhovz_i_in = rho_i_in * vz_in
                 p_i_in = p_in - p_e_in
                 e_i_in = p_i_in / (gasGamma - 1.) + 0.5 * (rhovx_i_in^2 + rhovy_i_in^2 + rhovz_i_in^2) / rho_i_in
               Pxx_i_in = p_i_in + vx_in^2 * rho_i_in
               Pyy_i_in = p_i_in + vy_in^2 * rho_i_in
               Pzz_i_in = p_i_in + vz_in^2 * rho_i_in
               Pxy_i_in = vx_in * vy_in * rho_i_in
               Pxz_i_in = vx_in * vz_in * rho_i_in
               Pyz_i_in = vy_in * vz_in * rho_i_in

                   B_in = math.sqrt(Bx_in^2 + By_in^2 + Bz_in^2)
                  cs_in = math.sqrt(gasGamma*p_in/rho_in)
                  vA_in = B_in/math.sqrt(mu0*rho_in)
                pmag_in = B_in^2/2/mu0
                beta_in = p_in/pmag_in
                  Ex_in = - vy_in*Bz_in + vz_in*By_in
                  Ey_in = - vz_in*Bx_in + vx_in*Bz_in
                  Ez_in = - vx_in*By_in + vy_in*Bx_in
                   dMin = math.min( (xup-xlo)/nx, (yup-ylo)/ny, (zup-zlo)/nz )

log("%30s = %g", "lightSpeed [km/s]", lightSpeed/1e3)
log("%30s = %g, %g, %g", "Dipole strength [nT]", Dx/Re^3/1e-9, Dy/Re^3/1e-9, Dz/Re^3/1e-9)
log("%30s = %g", "rho_in [amu/cm^3]", rho_in/Constants.PROTON_MASS/1e6)
log("%30s = %g, %g, %g", "v_in [km/s]", vx_in/1e3, vy_in/1e3, vz_in/1e3)
log("%30s = %g", "p_in [nPa]", p_in/1e-9)
log("%30s = %g, %g, %g", "B_in [nT]", Bx_in/1e-9, By_in/1e-9, Bz_in/1e-9)
log("%30s = %g", "cs_in [km/s]", cs_in/1e3)
log("%30s = %g", "vA_in [km/s]", vA_in/1e3)
log("%30s = %g", "pmag_in [nPa]", pmag_in/1e-9)
log("%30s = %g", "beta_in", beta_in)
log("%30s = %g", "r_ramp1 [Re]", r_ramp1/Re)
log("%30s = %g", "r_ramp2 [Re]", r_ramp2/Re)
log("%30s = %g", "stretch", stretch)
log("%30s = %g", "xmir [Re]", xmir/Re)
log("%30s = %g", "d_i_in [Re]", d_i_in/Re)
log("%30s = %g", "mass ratio", mass_ratio)
log("%30s = %g", "pressure ratio", pressure_ratio)
log("%30s = %g", "elcMass [me]", elcMass/Constants.ELECTRON_MASS)
log("%30s = %g", "ionMass [mp]", ionMass/Constants.PROTON_MASS)
log("%30s = %g", "elcCharge [e]", elcCharge/Constants.ELEMENTARY_CHARGE)
log("%30s = %g", "ionCharge [e]", ionCharge/Constants.ELEMENTARY_CHARGE)
log("%30s = %g", "ionMass/ionCharge [mp/e]", ionMass/ionCharge / (Constants.PROTON_MASS/Constants.ELEMENTARY_CHARGE))
log("%30s = %g", "elcMass/elcCharge [me/e]", elcMass/elcCharge / (Constants.ELECTRON_MASS/Constants.ELEMENTARY_CHARGE))
log("%30s = %g, %g, %g", "lower [Re]", xlo/Re, ylo/Re, zlo/Re)
log("%30s = %g, %g, %g", "upper [Re]", xup/Re, yup/Re, zup/Re)
log("%30s = %g, %g, %g", "cells", nx, ny, nz)

-----------------------
-- INITIAL CONDITION --
-----------------------
function calcRho(x,y,z)
   local rho = rho_in
   return rho
end

function calcP(x,y,z)
   local p = p_in
   return p
end

function calcV(x,y,z)
   local xx = x > 0 and x/stretch or x
   local r = math.sqrt(xx^2 + y^2 + z^2)
   local s = (r-r_ramp1)/(r_ramp2-r_ramp1)
   s = math.max(s, 0)
   s = math.min(s, 1)
   local vx = vx_in * s
   local vy = vy_in * s
   local vz = vz_in * s
   return vx, vy, vz
end

function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz, cut)
   if cut and x < xmir then
      return 0, 0, 0
   end
   local xx = (x - x0)
   local yy = (y - y0)
   local zz = (z - z0)
   local rr = math.sqrt(xx^2+yy^2+zz^2)
   local Bx = (3*xx*Dx*xx + 3*xx*Dy*yy + 3*xx*Dz*zz - Dx*rr^2) / rr^5
   local By = (3*yy*Dx*xx + 3*yy*Dy*yy + 3*yy*Dz*zz - Dy*rr^2) / rr^5
   local Bz = (3*zz*Dx*xx + 3*zz*Dy*yy + 3*zz*Dz*zz - Dz*rr^2) / rr^5
   return Bx, By, Bz
end

function staticB(x, y, z)
   local r2 = x^2+y^2+z^2
   if (r2 < r1^2) then
      return 0, 0, 0
   end
   local Bxs, Bys, Bzs = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, false)
   Bxs, Bys, Bzs = Bxs+Bxd0, Bys+Byd0, Bzs+Bzd0
   -- Bxs, Bys, Bzs = Bxs+Bx_in, Bys+By_in, Bzs+Bz_in
   return Bxs, Bys, Bzs
end

function totalB(x,y,z)
   local r2 = x^2+y^2+z^2
   if (r2 < r1^2) then
      return 0, 0, 0
   end
   local Bxt, Byt, Bzt = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, true)
   Bxt, Byt, Bzt = Bxt+Bxd0, Byt+Byd0, Bzt+Bzd0
   local Bxdm, Bydm, Bzdm = dipoleB(x, y, z, 2*xmir, 0, 0, -Dx, Dy, Dz, true)
   Bxt, Byt, Bzt = Bxt+Bxdm, Byt+Bydm, Bzt+Bzdm
   Bxt, Byt, Bzt = Bxt+Bx_in, Byt+By_in, Bzt+Bz_in
   return Bxt, Byt, Bzt
end

function init(x,y,z)
   local rho = calcRho(x,y,z)
   local vx,vy,vz = calcV(x,y,z)
   local p = calcP(x,y,z)

   local rho_e = rho / (1 + mass_ratio)
   local rho_i = rho - rho_e
   local rhovx_e = rho_e * vx
   local rhovy_e = rho_e * vy
   local rhovz_e = rho_e * vz
   local rhovx_i = rho_i * vx
   local rhovy_i = rho_i * vy
   local rhovz_i = rho_i * vz
   local p_e = p / (1 + pressure_ratio)
   local p_i = p - p_e

   local e_e = p_e / (gasGamma - 1.) + 0.5 * (rhovx_e^2 + rhovy_e^2 + rhovz_e^2) / rho_e
   local e_i = p_i / (gasGamma - 1.) + 0.5 * (rhovx_i^2 + rhovy_i^2 + rhovz_i^2) / rho_i

   local Pxx_e = p_e + rhovx_e^2 / rho_e
   local Pyy_e = p_e + rhovy_e^2 / rho_e
   local Pzz_e = p_e + rhovz_e^2 / rho_e
   local Pxy_e = rhovx_e * rhovy_e / rho_e
   local Pxz_e = rhovx_e * rhovz_e / rho_e
   local Pyz_e = rhovy_e * rhovz_e / rho_e

   local Pxx_i = p_i + rhovx_i^2 / rho_i
   local Pyy_i = p_i + rhovy_i^2 / rho_i
   local Pzz_i = p_i + rhovz_i^2 / rho_i
   local Pxy_i = rhovx_i * rhovy_i / rho_i
   local Pxz_i = rhovx_i * rhovz_i / rho_i
   local Pyz_i = rhovy_i * rhovz_i / rho_i

   local Bxt, Byt, Bzt = totalB(x, y, z)
   local Bxs, Bys, Bzs = staticB(x, y, z)
   local Bx, By, Bz = Bxt-Bxs, Byt-Bys, Bzt-Bzs

   local Ex = - vy*Bzt + vz*Byt
   local Ey = - vz*Bxt + vx*Bzt
   local Ez = - vx*Byt + vy*Bxt

   return rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
          rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
          Ex, Ey, Ez, Bx, By, Bz, 0, 0
end

function setStaticField(x, y, z)
   local Exs, Eys, Ezs = 0, 0, 0
   local Bxs, Bys, Bzs = staticB(x, y, z)
   return Exs, Eys, Ezs, Bxs, Bys, Bzs, 0, 0
end

function setInOutField(x, y, z)
   if (x^2 + y^2 + z^2 < r0^2) then
      return -1
   else
      return 1
   end
end


----------
-- GRID --
----------

local coordinateMap = nil
if useNonUniformGrid then

-- distribution of grid size centered at x0
   function f_dx5(x, x0, w, r, uniform_width)
   --[[
   Args:
      x    : relative coordinate in the range of [0,1]
      xx0  : center of dx distribution on [0,1] that is transformed from 
             physical domain [xlo, xup]
      width: width of Gaussian distribution
   Returns:
   --]]
      if (math.abs(x-x0) <= uniform_width) then
         return 1
      elseif (x-x0 > uniform_width) then
         local x1 = x0 + uniform_width
         return (1/r - (1/r-1) * math.exp( -.5 * ((x-x1)/w)^2 ))
      elseif (x-x0 < -uniform_width) then
         local x1 = x0 - uniform_width
         return (1/r - (1/r-1) * math.exp( -.5 * ((x-x1)/w)^2 ))
      end
   end

   function gentle_peak1(x, x1, w1, r1)
      return (1/r1 - (1/r1-1) * math.exp( -0.5 * ((x-x1)/w1)^2 ))
   end

   function gentle_peak2(x, x1, wl, rl, wr, rr)
      local w1 = (x>x1) and wr or wl
      local r1 = (x>x1) and rr or rl
      return (1/r1 - (1/r1-1) * math.exp( -0.5 * ((x-x1)/w1)^2 ))
   end

   function calcGrid(f_dx, xlo, xup, nx)
   --[[
   Args:
      xlo, xup: lower and upper bounds of physical domain
      nx   : number of real cells
   Returns:
      xx_nc: node-center compuational domain coordinates; 'ghost' points might
             go out of [0, 1]
      x_nc : node-center physical domain coordinates
   Both xx_nc and x_c have nx+1+2*sw points. xx_nc is uniform and 1-on-1 maps
   to x_nc.
   --]]
      local sw = 2
      local sw1 = sw + 100 -- wider stencile for xx_nc to include more possible xx values
      local xx_nc = {}
      xx_nc[sw1+1] = 0
      for j=0,nx+sw1-1 do
         local i = j+sw1+1
         xx_nc[i+1] = xx_nc[i] + f_dx((j+.5)/nx)
      end
      for j=0,-sw1+1,-1 do
         local i = j+sw1+1
         xx_nc[i-1] = xx_nc[i] - f_dx((j-.5)/nx)
      end
      fac = 1/xx_nc[nx+sw1+1]
      for j=-sw1,nx+sw1 do
         local i = j+sw1+1
         xx_nc[i] = xx_nc[i]*fac
      end
      local x_nc = {}
      for j=-sw1,nx+sw1 do
         local i = j+sw1+1
         x_nc[i] = xlo + (xup-xlo) * j/nx
      end
      return xx_nc, x_nc
   end
   function linearRefine(xx, xx_nc, x_nc, display)
      local num = table.getn(xx_nc)
      local x
      if (xx <= xx_nc[1]) then
         x = x_nc[1]
      elseif (xx >= xx_nc[num]) then
         x = x_nc[num]
      end
      for i = 2,num do
         if (xx == xx_nc[i]) then
            x = x_nc[i]
            break
         elseif (xx < xx_nc[i]) then -- xx is between xx_nc[i-1] and xx_nc[i]
            x = x_nc[i-1] + (x_nc[i]-x_nc[i-1])*(xx-xx_nc[i-1])/(xx_nc[i]-xx_nc[i-1])
            break
         end
      end
      if display then
         log("%g %g", xx, x)
      end
      return x
   end

   --[[
   function f_dx(x)
      return gentle_peak2(x, (x1-xlo)/(xup-xlo), wx1l/(xup-xlo), rx1l, wx1r/(xup-xlo), rx1r) 
           + gentle_peak2(x, (x2-xlo)/(xup-xlo), wx2l/(xup-xlo), rx2l, wx2r/(xup-xlo), rx2r)
   end
   --]]
   x0 = 0
   function f_dx(x)
      return f_dx5(x, (x0 - xlo) / (xup - xlo), wx/(xup-xlo),
         rx, fwx/(xup-xlo))
   end
   xx_nc, x_nc = calcGrid(f_dx, xlo, xup, nx)

   y0 = 0 -- minimum dy location
   function f_dy(y)
      return f_dx5(y, (y0 - ylo) / (yup - ylo), wy/(yup-ylo),
         ry, fwy/(yup-ylo))
   end
   yy_nc, y_nc = calcGrid(f_dy, ylo, yup, ny)

   z0 = 0 -- minimum dz location
   function f_dz(z)
      return f_dx5(z, (z0 - zlo) / (zup - zlo), wz/(zup-zlo),
         rz, fwz/(zup-zlo))
   end
   zz_nc, z_nc = calcGrid(f_dz, zlo, zup, nz)

   coordinateMap = {
      function (xx)
         return linearRefine(xx, xx_nc, x_nc)
      end,                                  
      function (yy)                         
         return linearRefine(yy, yy_nc, y_nc)
      end,                                  
      function (zz)                         
         return linearRefine(zz, zz_nc, z_nc)
      end,
   }

   lower = {0, 0, 0}
   upper = {1, 1, 1}
end


local function build(self, tbl)
   local tEnd = tonumber(tbl.tEnd)

   local momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = lower,
   upper = upper,
   cells = cells,
   coordinateMap = coordinateMap,
   timeStepper = "fvDimSplit",

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      evolve = true,

      -- initial conditions
      init = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]

         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
               rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
               Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x, y, z)

         return rho_e, rhovx_e, rhovy_e, rhovz_e, e_e
      end,

      bcx = { Euler.bcConst(rho_e_in, rhovx_e_in, rhovy_e_in, rhovz_e_in, e_e_in),
              Euler.bcCopy },
      bcy = { Euler.bcCopy, Euler.bcCopy },
      bcz = { Euler.bcCopy, Euler.bcCopy },

      hasSsBnd = true,
      inOutFunc = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = { Moments.Species.bcCopy },
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,

      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      evolve = true,

      -- initial conditions
      init = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
   
         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
               rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
               Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x, y, z)
   
         return rho_i, rhovx_i, rhovy_i, rhovz_i, e_i

      end,

      bcx = { Euler.bcConst(rho_i_in, rhovx_i_in, rhovy_i_in, rhovz_i_in, e_i_in),
              Euler.bcCopy },
      bcy = { Euler.bcCopy, Euler.bcCopy },
      bcz = { Euler.bcCopy, Euler.bcCopy },

      hasSsBnd = true,
      inOutFunc = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = { Moments.Species.bcCopy },
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      evolve = true,

      init = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]

         local rho_e, rhovx_e, rhovy_e, rhovz_e, e_e,
               rho_i, rhovx_i, rhovy_i, rhovz_i, e_i,
               Ex, Ey, Ez, Bx, By, Bz, phiE, phiB = init(x, y, z)

         return Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
      end,

      bcx = { Moments.Field.bcConst(Ex_in, Ey_in, Ez_in, Bx_in, By_in, Bz_in, 0, 0),
              Moments.Field.bcCopy },
      bcy = { Moments.Field.bcCopy, Moments.Field.bcCopy },
      bcz = { Moments.Field.bcCopy, Moments.Field.bcCopy },

      hasSsBnd = true,
      inOutFunc = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = { Moments.Species.bcReflect },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
      hasStaticField = true,
      staticEmFunction = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setStaticField(x, y, z)
      end
   },

}
return momentApp
end


local App = Proto()

function App:init(tbl)
end

function App:fullInit(tbl)
   self.app = build(self, tbl)
end

return App

