--------------
-- PREAMBLE --
--------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local TenMoment = require "Eq.TenMoment"
local Consts = require "Constants"
local Logger = require "Logger"

local logger = Logger {logToFile = True}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

local sqrt = math.sqrt
local min = math.min
local max = math.max
local abs = math.abs
local exp = math.exp
local floor = math.floor

----------------
-- PARAMETERS --
-- SI units   --
----------------
local gasGamma = 5 / 3
local lightSpeed = Consts.SPEED_OF_LIGHT / 100
local mu0 = Consts.MU0
local epsilon0 = 1 / mu0 / (lightSpeed ^ 2)

local R0 = 2634.1e3 -- Planet radius
-- Planet dipole moment
local Dx = -18e-9 * R0 ^ 3
local Dy = 51.8e-9 * R0 ^ 3
local Dz = -716.8e-9 * R0 ^ 3
local r1 = 0.5 * R0 -- Zero B field when r < r1
local r0 = R0 -- Radius of the inner boundary

-- Inflows condition
local rho_in = 56e6 * Consts.PROTON_MASS
local vx_in = 140e3
local vy_in = 0
local vz_in = 0
local p_in = 3.8e-9
local Bx_in = 0
local By_in = -6e-9
local Bz_in = -77e-9

-- Mirror-dipole setup
local r_ramp1 = 2 * R0
local r_ramp2 = 2.5 * R0
local stretch = 1
local xmir = -2.578125 * R0

-- Kinetic parameters
local mi = 14 * Consts.PROTON_MASS
local di_in = 0.2 * R0
local mi__me = 25
local pi__pe = 5
local me = mi / mi__me
local qi = mi / di_in / sqrt(mu0 * rho_in)
local qe = -qi
local charge = {qe, qi}
local mass = {me, mi}
local de_in = di_in * sqrt(me / mi)

-- Domain and grid
local xlo, xup, nx = -24 * R0, 24 * R0, 96
local ylo, yup, ny = -24 * R0, 24 * R0, 96
local zlo, zup, nz = -24 * R0, 24 * R0, 96
local lower = {xlo, ylo, zlo}
local upper = {xup, yup, zup}
local cells = {nx, ny, nz}
-- Nonuniform grid where grid size is a Flattop Gaussian function of location.
-- w: Gaussian region width.
-- fw: Flat region width.
-- r: Grid ratio, grid_max/grid_min.
local useNonUniformGrid = true
local wx, fwx, rx = 1.9 * R0, 2.25 * R0, 50
local wy, fwy, ry = 2. * R0, 2. * R0, 50
local wz, fwz, rz = 2. * R0, 2. * R0, 50

-- Computational constants.
local cfl = 0.9
local cflm = 1.1 * cfl
local limiter = "monotonized-centered"

-- I/O control
local tEnd = 2700
local tFrame = 60
local nFrame = floor(tEnd / tFrame)

-- Derived parameters
local rhoe_in = rho_in / (1 + mi__me)
local rhovxe_in = rhoe_in * vx_in
local rhovye_in = rhoe_in * vy_in
local rhovze_in = rhoe_in * vz_in
local pe_in = p_in / (1 + pi__pe)
local Pxxe_in = pe_in + vx_in ^ 2 * rhoe_in
local Pyye_in = pe_in + vy_in ^ 2 * rhoe_in
local Pzze_in = pe_in + vz_in ^ 2 * rhoe_in
local Pxye_in = vx_in * vy_in * rhoe_in
local Pxze_in = vx_in * vz_in * rhoe_in
local Pyze_in = vy_in * vz_in * rhoe_in

local rhoi_in = rho_in - rhoe_in
local rhovxi_in = rhoi_in * vx_in
local rhovyi_in = rhoi_in * vy_in
local rhovzi_in = rhoi_in * vz_in
local pi_in = p_in - pe_in
local Pxxi_in = pi_in + vx_in ^ 2 * rhoi_in
local Pyyi_in = pi_in + vy_in ^ 2 * rhoi_in
local Pzzi_in = pi_in + vz_in ^ 2 * rhoi_in
local Pxyi_in = vx_in * vy_in * rhoi_in
local Pxzi_in = vx_in * vz_in * rhoi_in
local Pyzi_in = vy_in * vz_in * rhoi_in

local B_in = sqrt(Bx_in ^ 2 + By_in ^ 2 + Bz_in ^ 2)
local cs_in = sqrt(gasGamma * p_in / rho_in)
local vA_in = B_in / sqrt(mu0 * rho_in)
local pmag_in = B_in ^ 2 / 2 / mu0
local beta_in = p_in / pmag_in
local Ex_in = -vy_in * Bz_in + vz_in * By_in
local Ey_in = -vz_in * Bx_in + vx_in * Bz_in
local Ez_in = -vx_in * By_in + vy_in * Bx_in

log("%30s = %g", "lightSpeed [km/s]", lightSpeed / 1e3)
log("%30s = %g, %g, %g", "Dipole strength [nT]", Dx / R0 ^ 3 / 1e-9,
    Dy / R0 ^ 3 / 1e-9, Dz / R0 ^ 3 / 1e-9)
log("%30s = %g", "rho_in [amu/cm^3]", rho_in / Consts.PROTON_MASS / 1e6)
log("%30s = %g, %g, %g", "v_in [km/s]", vx_in / 1e3, vy_in / 1e3, vz_in / 1e3)
log("%30s = %g", "p_in [nPa]", p_in / 1e-9)
log("%30s = %g, %g, %g", "B_in [nT]", Bx_in / 1e-9, By_in / 1e-9, Bz_in / 1e-9)
log("%30s = %g", "cs_in [km/s]", cs_in / 1e3)
log("%30s = %g", "vA_in [km/s]", vA_in / 1e3)
log("%30s = %g", "pmag_in [nPa]", pmag_in / 1e-9)
log("%30s = %g", "beta_in", beta_in)
log("%30s = %g", "mirror-dipole r_ramp1 [R0]", r_ramp1 / R0)
log("%30s = %g", "mirror-dipole r_ramp2 [R0]", r_ramp2 / R0)
log("%30s = %g", "mirror-dipole stretch     ", stretch)
log("%30s = %g", "mirror-dipole xmir    [R0]", xmir / R0)
log("%30s = %g", "di_in [R0]", di_in / R0)
log("%30s = %g", "ion-electron mass ratio", mi__me)
log("%30s = %g", "ion-electron pressure ratio", pi__pe)
log("%30s = %g", "me [true me]", me / Consts.ELECTRON_MASS)
log("%30s = %g", "mi [true mp]", mi / Consts.PROTON_MASS)
log("%30s = %g", "qe [true e]", qe / Consts.ELEMENTARY_CHARGE)
log("%30s = %g", "qi [true e]", qi / Consts.ELEMENTARY_CHARGE)
log("%30s = %g", "mi/qi [true mp/e]",
    mi / qi / (Consts.PROTON_MASS / Consts.ELEMENTARY_CHARGE))
log("%30s = %g", "me/qe [true me/e]",
    me / qe / (Consts.ELECTRON_MASS / Consts.ELEMENTARY_CHARGE))
log("%30s = %g, %g, %g", "lower [R0]", xlo / R0, ylo / R0, zlo / R0)
log("%30s = %g, %g, %g", "upper [R0]", xup / R0, yup / R0, zup / R0)
log("%30s = %g, %g, %g", "cells", nx, ny, nz)

-----------------------
-- INITIAL CONDITION --
-----------------------
local function calcRho(x, y, z)
   local rho = rho_in
   return rho
end

local function calcP(x, y, z)
   local p = p_in
   return p
end

local function calcV(x, y, z)
   local xx = x > 0 and x / stretch or x
   local r = sqrt(xx ^ 2 + y ^ 2 + z ^ 2)
   local s = (r - r_ramp1) / (r_ramp2 - r_ramp1)
   s = max(s, 0)
   s = min(s, 1)
   local vx = vx_in * s
   local vy = vy_in * s
   local vz = vz_in * s
   return vx, vy, vz
end

local function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz, cut)
   if cut and x < xmir then return 0, 0, 0 end
   local xx = (x - x0)
   local yy = (y - y0)
   local zz = (z - z0)
   local rr = sqrt(xx ^ 2 + yy ^ 2 + zz ^ 2)
   local Bx =
       (3 * xx * Dx * xx + 3 * xx * Dy * yy + 3 * xx * Dz * zz - Dx * rr ^ 2) /
           rr ^ 5
   local By =
       (3 * yy * Dx * xx + 3 * yy * Dy * yy + 3 * yy * Dz * zz - Dy * rr ^ 2) /
           rr ^ 5
   local Bz =
       (3 * zz * Dx * xx + 3 * zz * Dy * yy + 3 * zz * Dz * zz - Dz * rr ^ 2) /
           rr ^ 5
   return Bx, By, Bz
end

local function staticB(x, y, z)
   local r2 = x ^ 2 + y ^ 2 + z ^ 2
   if (r2 < r1 ^ 2) then return 0, 0, 0 end
   local Bxs, Bys, Bzs = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, false)
   Bxs, Bys, Bzs = Bxs + Bxd0, Bys + Byd0, Bzs + Bzd0
   -- Bxs, Bys, Bzs = Bxs+Bx_in, Bys+By_in, Bzs+Bz_in
   return Bxs, Bys, Bzs
end

local function totalB(x, y, z)
   local r2 = x ^ 2 + y ^ 2 + z ^ 2
   if (r2 < r1 ^ 2) then return 0, 0, 0 end
   local Bxt, Byt, Bzt = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, true)
   Bxt, Byt, Bzt = Bxt + Bxd0, Byt + Byd0, Bzt + Bzd0
   local Bxdm, Bydm, Bzdm = dipoleB(x, y, z, 2 * xmir, 0, 0, -Dx, Dy, Dz, true)
   Bxt, Byt, Bzt = Bxt + Bxdm, Byt + Bydm, Bzt + Bzdm
   Bxt, Byt, Bzt = Bxt + Bx_in, Byt + By_in, Bzt + Bz_in
   return Bxt, Byt, Bzt
end

local function init(xn, name)
   local x, y, z = xn[1], xn[2], xn[3]

   local rho = calcRho(x, y, z)
   local vx, vy, vz = calcV(x, y, z)
   local p = calcP(x, y, z)

   local rhoe = rho / (1 + mi__me)
   local rhoi = rho - rhoe
   local rhovxe = rhoe * vx
   local rhovye = rhoe * vy
   local rhovze = rhoe * vz
   local rhovxi = rhoi * vx
   local rhovyi = rhoi * vy
   local rhovzi = rhoi * vz
   local pe = p / (1 + pi__pe)
   local pi = p - pe

   local Pxxe = pe + rhovxe ^ 2 / rhoe
   local Pyye = pe + rhovye ^ 2 / rhoe
   local Pzze = pe + rhovze ^ 2 / rhoe
   local Pxye = rhovxe * rhovye / rhoe
   local Pxze = rhovxe * rhovze / rhoe
   local Pyze = rhovye * rhovze / rhoe

   local Pxxi = pi + rhovxi ^ 2 / rhoi
   local Pyyi = pi + rhovyi ^ 2 / rhoi
   local Pzzi = pi + rhovzi ^ 2 / rhoi
   local Pxyi = rhovxi * rhovyi / rhoi
   local Pxzi = rhovxi * rhovzi / rhoi
   local Pyzi = rhovyi * rhovzi / rhoi

   local Bxt, Byt, Bzt = totalB(x, y, z)
   local Bxs, Bys, Bzs = staticB(x, y, z)
   local Bx, By, Bz = Bxt - Bxs, Byt - Bys, Bzt - Bzs

   local Ex = -vy * Bzt + vz * Byt
   local Ey = -vz * Bxt + vx * Bzt
   local Ez = -vx * Byt + vy * Bxt

   if name == "elc" then
      return rhoe, rhovxe, rhovye, rhovze, Pxxe, Pxye, Pxze, Pyye, Pyze, Pzze
   elseif name == "ion" then
      return rhoi, rhovxi, rhovyi, rhovzi, Pxxi, Pxyi, Pxzi, Pyyi, Pyzi, Pzzi
   elseif name == "field" then
      return Ex, Ey, Ez, Bx, By, Bz, 0, 0
   else
      assert(0)
   end
end

local function setStaticField(x, y, z)
   local Exs, Eys, Ezs = 0, 0, 0
   local Bxs, Bys, Bzs = staticB(x, y, z)
   return Exs, Eys, Ezs, Bxs, Bys, Bzs, 0, 0
end

local function setInOutField(x, y, z)
   if (x ^ 2 + y ^ 2 + z ^ 2 < r0 ^ 2) then
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
      if (abs(x - x0) <= uniform_width) then
         return 1
      elseif (x - x0 > uniform_width) then
         local x1 = x0 + uniform_width
         return (1 / r - (1 / r - 1) * exp(-.5 * ((x - x1) / w) ^ 2))
      elseif (x - x0 < -uniform_width) then
         local x1 = x0 - uniform_width
         return (1 / r - (1 / r - 1) * exp(-.5 * ((x - x1) / w) ^ 2))
      end
   end

   function gentle_peak1(x, x1, w1, r1)
      return (1 / r1 - (1 / r1 - 1) * exp(-0.5 * ((x - x1) / w1) ^ 2))
   end

   function gentle_peak2(x, x1, wl, rl, wr, rr)
      local w1 = (x > x1) and wr or wl
      local r1 = (x > x1) and rr or rl
      return (1 / r1 - (1 / r1 - 1) * exp(-0.5 * ((x - x1) / w1) ^ 2))
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
      xx_nc[sw1 + 1] = 0
      for j = 0, nx + sw1 - 1 do
         local i = j + sw1 + 1
         xx_nc[i + 1] = xx_nc[i] + f_dx((j + .5) / nx)
      end
      for j = 0, -sw1 + 1, -1 do
         local i = j + sw1 + 1
         xx_nc[i - 1] = xx_nc[i] - f_dx((j - .5) / nx)
      end
      fac = 1 / xx_nc[nx + sw1 + 1]
      for j = -sw1, nx + sw1 do
         local i = j + sw1 + 1
         xx_nc[i] = xx_nc[i] * fac
      end
      local x_nc = {}
      for j = -sw1, nx + sw1 do
         local i = j + sw1 + 1
         x_nc[i] = xlo + (xup - xlo) * j / nx
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
      for i = 2, num do
         if (xx == xx_nc[i]) then
            x = x_nc[i]
            break
         elseif (xx < xx_nc[i]) then -- xx is between xx_nc[i-1] and xx_nc[i]
            x = x_nc[i - 1] + (x_nc[i] - x_nc[i - 1]) * (xx - xx_nc[i - 1]) /
                    (xx_nc[i] - xx_nc[i - 1])
            break
         end
      end
      if display then log("%g %g", xx, x) end
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
      return f_dx5(x, (x0 - xlo) / (xup - xlo), wx / (xup - xlo), rx,
                   fwx / (xup - xlo))
   end
   xx_nc, x_nc = calcGrid(f_dx, xlo, xup, nx)

   y0 = 0 -- minimum dy location
   function f_dy(y)
      return f_dx5(y, (y0 - ylo) / (yup - ylo), wy / (yup - ylo), ry,
                   fwy / (yup - ylo))
   end
   yy_nc, y_nc = calcGrid(f_dy, ylo, yup, ny)

   z0 = 0 -- minimum dz location
   function f_dz(z)
      return f_dx5(z, (z0 - zlo) / (zup - zlo), wz / (zup - zlo), rz,
                   fwz / (zup - zlo))
   end
   zz_nc, z_nc = calcGrid(f_dz, zlo, zup, nz)

   coordinateMap = {
      function(xx) return linearRefine(xx, xx_nc, x_nc) end,
      function(yy) return linearRefine(yy, yy_nc, y_nc) end,
      function(zz) return linearRefine(zz, zz_nc, z_nc) end
   }

   lower = {0, 0, 0}
   upper = {1, 1, 1}
end

---------
-- APP --
---------
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
      charge = qe,
      mass = me,
      equation = TenMoment {},
      equationInv = TenMoment {numericalFlux = "lax"},
      init = function(t, xn) return init(xn, "elc") end,
      bcx = {
         TenMoment.bcConst(rhoe_in, rhovxe_in, rhovye_in, rhovze_in, Pxxe_in,
                           Pxye_in, Pxze_in, Pyye_in, Pyze_in, Pzze_in),
         Moments.Species.bcCopy
      },
      bcy = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      bcz = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      hasSsBnd = true,
      inOutFunc = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = {Moments.Species.bcCopy}
   },

   -- ions
   ion = Moments.Species {
      charge = qi,
      mass = mi,
      equation = TenMoment {gasGamma = gasGamma},
      equationInv = TenMoment {gasGamma = gasGamma, numericalFlux = "lax"},
      init = function(t, xn) return init(xn, "ion") end,
      bcx = {
         TenMoment.bcConst(rhoi_in, rhovxi_in, rhovyi_in, rhovzi_in, Pxxi_in,
                           Pxyi_in, Pxzi_in, Pyyi_in, Pyzi_in, Pzzi_in),
         Moments.Species.bcCopy
      },
      bcy = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      bcz = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      hasSsBnd = true,
      inOutFunc = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = {Moments.Species.bcCopy}
   },

   field = Moments.Field {
      epsilon0 = 1.0,
      mu0 = 1.0,
      init = function(t, xn) return init(xn, "field") end,
      bcx = {
         Moments.Field.bcConst(Ex_in, Ey_in, Ez_in, Bx_in, By_in, Bz_in, 0, 0),
         Moments.Field.bcCopy
      },
      bcy = {Moments.Field.bcCopy, Moments.Field.bcCopy},
      bcz = {Moments.Field.bcCopy, Moments.Field.bcCopy},
      inOutFunc = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setInOutField(x, y, z)
      end,
      ssBc = {Moments.Species.bcReflect}
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
      hasStaticField = true,
      staticEmFunction = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return setStaticField(x, y, z)
      end
   }

}

---------
-- RUN --
---------
momentApp:run()
