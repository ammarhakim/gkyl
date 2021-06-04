-----------------------------------------------------------------
-- Two-fluid modeling of mangetosphere-solar wind interaction. --
-----------------------------------------------------------------
--------------
-- PREAMBLE --
--------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local TenMoment = require "Eq.TenMoment"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Consts = require "Constants"
local Logger = require "Logger"

local logger = Logger {logToFile = True}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

----------------
-- PARAMETERS --
-- SI units   --
----------------
local gasGamma = 5 / 3
local lightSpeed = Consts.SPEED_OF_LIGHT / 100
local mu0 = Consts.MU0
local epsilon0 = 1 / mu0 / (lightSpeed ^ 2)

local R0 = 2634.1e3 -- Planet radius.
-- Planet dipole moment.
local Dx = -18e-9 * R0 ^ 3
local Dy = 51.8e-9 * R0 ^ 3
local Dz = -716.8e-9 * R0 ^ 3
local r1 = 0.5 * R0 -- B field will be set to zero where r < r1.
local rInOut = R0 -- Radius of the inner boundary.

-- Inflow condition.
local rho_in = 56e6 * Consts.PROTON_MASS
local vx_in = 140e3
local vy_in = 0
local vz_in = 0
local p_in = 3.8e-9
local Bx_in = 0
local By_in = -6e-9
local Bz_in = -77e-9

-- Mirror-dipole (mirdip) setup.
local mirdip_rRamp1 = 2 * R0
local mirdip_rRamp2 = 2.5 * R0
local mirdip_stretch = 1
local mirdip_xMirror = -2.578125 * R0

-- Kinetic parameters.
local mi = 14 * Consts.PROTON_MASS
local di_in = 0.2 * R0
local mi__me = 25
local pi__pe = 5
local me = mi / mi__me
local qi = mi / di_in / math.sqrt(mu0 * rho_in)
local qe = -qi
local charge = {qe, qi}
local mass = {me, mi}
local de_in = di_in * math.sqrt(me / mi)

-- Domain.
local xlo, xup, Nx = -24 * R0, 24 * R0, 96
local ylo, yup, Ny = -24 * R0, 24 * R0, 96
local zlo, zup, Nz = -24 * R0, 24 * R0, 96
-- Grid.
local lower = {xlo, ylo, zlo}
local upper = {xup, yup, zup}
local cells = {Nx, Ny, Nz}
local decompCuts = nil
local useNonUniformGrid = true
-- Using two math.tanh functions to rampd down and up the grid sizes (U-shape).
-- xl, wxl: Floor and transition-layer width of the left math.tanh (ramping down dx.)
-- xr, wxr: Floor and transition-layer width of the right math.tanh (ramping up dx.)
-- sx: Shift all grids up/down to avoid tiny grids and to somewhat change the
-- math.max(dx)/math.min(dx) ratio.
local xl, wxl, xr, wxr, sx = -5 * R0, 3 * R0, 5 * R0, 3 * R0, 0.01
local yl, wyl, yr, wyr, sy = -5 * R0, 3 * R0, 5 * R0, 3 * R0, 0.01
local zl, wzl, zr, wzr, sz = -5 * R0, 3 * R0, 5 * R0, 3 * R0, 0.01

-- Computational constants.
local cfl = 0.9
local cflm = 1.1 * cfl
local limiter = "monotonized-centered"

-- I/O control.
local tEnd = 2700
local tFrame = 60
local nFrame = math.floor(tEnd / tFrame)

-- Derived parameters.
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

local B_in = math.sqrt(Bx_in ^ 2 + By_in ^ 2 + Bz_in ^ 2)
local cs_in = math.sqrt(gasGamma * p_in / rho_in)
local vA_in = B_in / math.sqrt(mu0 * rho_in)
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
log("%30s = %g", "mirdip_rRamp1  [R0]", mirdip_rRamp1 / R0)
log("%30s = %g", "mirdip_rRamp2  [R0]", mirdip_rRamp2 / R0)
log("%30s = %g", "mirdip_stretch     ", mirdip_stretch)
log("%30s = %g", "mirdip_xMirror [R0]", mirdip_xMirror / R0)
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
log("%30s = %g, %g, %g", "cells", Nx, Ny, Nz)

-----------------------
-- INITIAL CONDITION --
-----------------------
local function initRho(x, y, z)
   local rho = rho_in
   return rho
end

local function initP(x, y, z)
   local p = p_in
   return p
end

local function initV(x, y, z)
   local xx = x > 0 and x / mirdip_stretch or x
   local r = math.sqrt(xx ^ 2 + y ^ 2 + z ^ 2)
   local s = (r - mirdip_rRamp1) / (mirdip_rRamp2 - mirdip_rRamp1)
   s = math.max(s, 0)
   s = math.min(s, 1)
   local vx = vx_in * s
   local vy = vy_in * s
   local vz = vz_in * s
   return vx, vy, vz
end

local function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz, cut)
   if cut and x < mirdip_xMirror then return 0, 0, 0 end
   local xx, yy, zz = x - x0, y - y0, z - z0
   local rr = math.sqrt(xx ^ 2 + yy ^ 2 + zz ^ 2)
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
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, false)
   return Bxd0, Byd0, Bzd0
end

local function totalB(x, y, z)
   local r2 = x ^ 2 + y ^ 2 + z ^ 2
   if (r2 < r1 ^ 2) then return 0, 0, 0 end

   local Bxt, Byt, Bzt = 0, 0, 0

   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, true)
   Bxt, Byt, Bzt = Bxt + Bxd0, Byt + Byd0, Bzt + Bzd0

   local Bxdm, Bydm, Bzdm = dipoleB(x, y, z, 2 * mirdip_xMirror, 0, 0, -Dx, Dy,
                                    Dz, true)
   Bxt, Byt, Bzt = Bxt + Bxdm, Byt + Bydm, Bzt + Bzdm

   Bxt, Byt, Bzt = Bxt + Bx_in, Byt + By_in, Bzt + Bz_in
   return Bxt, Byt, Bzt
end

local rho_v_p_to_10m = function(q, rho, vx, vy, vz, p)
   q[1] = rho
   q[2] = rho * vx
   q[3] = rho * vy
   q[4] = rho * vz
   q[5] = rho * vx * vx + p
   q[6] = rho * vx * vy
   q[7] = rho * vx * vz
   q[8] = rho * vy * vy + p
   q[9] = rho * vy * vz
   q[10] = rho * vz * vz + p
end

local function initElc(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]

   local rhoe = initRho(x, y, z) / (1 + mi__me)
   local vx, vy, vz = initV(x, y, z)
   local pe = initP(x, y, z) / (1 + pi__pe)

   local q = {}
   rho_v_p_to_10m(q, rhoe, vx, vy, vz, pe)

   return unpack(q)
end

local function initIon(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]

   local rhoi = initRho(x, y, z) * mi__me / (1 + mi__me)
   local vx, vy, vz = initV(x, y, z)
   local pi = initP(x, y, z) * pi__pe / (1 + pi__pe)

   local q = {}
   rho_v_p_to_10m(q, rhoi, vx, vy, vz, pi)

   return unpack(q)
end

local function initField(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]

   local Bxt, Byt, Bzt = totalB(x, y, z)
   local Bxs, Bys, Bzs = staticB(x, y, z)
   local Bx, By, Bz = Bxt - Bxs, Byt - Bys, Bzt - Bzs

   local vx, vy, vz = initV(x, y, z)
   local Ex = -vy * Bzt + vz * Byt
   local Ey = -vz * Bxt + vx * Bzt
   local Ez = -vx * Byt + vy * Bxt

   return Ex, Ey, Ez, Bx, By, Bz, 0, 0
end

local function staticEmFunction(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local Exs, Eys, Ezs = 0, 0, 0
   local Bxs, Bys, Bzs = staticB(x, y, z)
   return Exs, Eys, Ezs, Bxs, Bys, Bzs, 0, 0
end

local function inOutFunc(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if (x ^ 2 + y ^ 2 + z ^ 2 < rInOut ^ 2) then
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
   -- User-defined functions that map computational-grid coordinates to
   -- computational-grid cell sizes. The physical coordinates will be computed
   -- accordingly.
   local Lx = xup - xlo
   local _xl, _wxl = (xl - xlo) / Lx, wxl / Lx
   local _xr, _wxr = (xr - xlo) / Lx, wxr / Lx
   local _x2dx = function(_x)
      return 2 - math.tanh((_x - _xl) / _wxl) + math.tanh((_x - _xr) / _wxr) + sx
   end

   local Ly = yup - ylo
   local _yl, _wyl = (yl - ylo) / Ly, wyl / Ly
   local _yr, _wyr = (yr - ylo) / Ly, wyr / Ly
   local _y2dy = function(_y)
      return 2 - math.tanh((_y - _yl) / _wyl) + math.tanh((_y - _yr) / _wyr) + sy
   end

   local Lz = zup - zlo
   local _zl, _wzl = (zl - zlo) / Lz, wzl / Lz
   local _zr, _wzr = (zr - zlo) / Lz, wzr / Lz
   local _z2dz = function(_z)
      return 2 - math.tanh((_z - _zl) / _wzl) + math.tanh((_z - _zr) / _wzr) + sz
   end

   -- Get a mapping function from computational node-center coordinates to
   -- physical node-center coordinates.
   -- _x2dx: maps cell index to grid size in the computational grid.
   local get_c2p = function(x0, Lx, Nx, _x2dx)
      -- Node-center computational coordinates.
      local xncComput = {0}
      for i = 2, Nx + 1 do
         xncComput[i] = xncComput[i - 1] + _x2dx((i - 0.5) / Nx)
      end

      -- Node-center physical coordinates.
      local xncPhys = {}
      for i = 1, Nx + 1 do
         xncPhys[i] = x0 + xncComput[i] * Lx / xncComput[Nx + 1]
      end

      xncComput = nil
      return function(xncComput)
         local idx = math.floor(xncComput * Nx + 1.5)
         idx = math.max(1, idx)
         idx = math.min(Nx + 1, idx)
         return xncPhys[idx]
      end
   end

   coordinateMap = {
      get_c2p(xlo, Lx, Nx, _x2dx), get_c2p(ylo, Ly, Ny, _y2dy),
      get_c2p(zlo, Lz, Nz, _z2dz)
   }
   lower = {0, 0, 0}
   upper = {1, 1, 1}
end

------------------------
-- BOUNDARY CONDITION --
------------------------
local bcInflow_elc = TenMoment.bcConst(rhoe_in, rhovxe_in, rhovye_in, rhovze_in,
                                       Pxxe_in, Pxye_in, Pxze_in, Pyye_in,
                                       Pyze_in, Pzze_in)
local bcInflow_ion = TenMoment.bcConst(rhoi_in, rhovxi_in, rhovyi_in, rhovzi_in,
                                       Pxxi_in, Pxyi_in, Pxzi_in, Pyyi_in,
                                       Pyzi_in, Pzzi_in)

local bcInflow_field = {
   function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
      qbc[1] = Ex_in
      qbc[2] = Ey_in
      qbc[3] = Ez_in

      local x, y, z = xcOut[1], xcOut[2], xcOut[3]
      local Bxs, Bys, Bzs = staticB(x, y, z)
      qbc[4] = Bx_in - Bxs
      qbc[5] = By_in - Bys
      qbc[6] = Bz_in - Bzs

      qbc[7] = qin[7]
      qbc[8] = qin[8]
   end
}

local bcInner_elc = {
   function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
      rho_v_p_to_10m(qbc, rhoe_in, 0, 0, 0, pe_in)
   end
}

local bcInner_ion = {
   function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
      rho_v_p_to_10m(qbc, rhoi_in, 0, 0, 0, pi_in)
   end
}

local bcInner_field = {
   function(dir, tm, idxIn, qin, qbc, xcOut, xcIn)
      qbc[1] = 0
      qbc[2] = 0
      qbc[3] = 0

      qbc[4] = qin[4]
      qbc[5] = qin[5]
      qbc[6] = qin[6]

      qbc[7] = qin[7]
      qbc[8] = qin[8]
   end
}

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
   decompCuts = decompCuts,
   timeStepper = "fvDimSplit",

   elc = Moments.Species {
      charge = qe,
      mass = me,
      equation = TenMoment {},
      equationInv = TenMoment {numericalFlux = "lax"},
      init = initElc,
      bcx = {bcInflow_elc, Moments.Species.bcCopy},
      bcy = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      bcz = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      hasSsBnd = true,
      inOutFunc = inOutFunc,
      ssBc = {bcInner_elc}
   },

   ion = Moments.Species {
      charge = qi,
      mass = mi,
      equation = TenMoment {gasGamma = gasGamma},
      equationInv = TenMoment {gasGamma = gasGamma, numericalFlux = "lax"},
      init = initIon,
      bcx = {bcInflow_ion, Moments.Species.bcCopy},
      bcy = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      bcz = {Moments.Species.bcCopy, Moments.Species.bcCopy},
      hasSsBnd = true,
      inOutFunc = inOutFunc,
      ssBc = {bcInner_ion}
   },

   field = Moments.Field {
      epsilon0 = epsilon0,
      mu0 = mu0,
      init = initField,
      bcx = {bcInflow_field, Moments.Field.bcCopy},
      bcy = {Moments.Field.bcCopy, Moments.Field.bcCopy},
      bcz = {Moments.Field.bcCopy, Moments.Field.bcCopy},
      hasSsBnd = true,
      inOutFunc = inOutFunc,
      ssBc = {bcInner_field}
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
      hasStaticField = true,
      staticEmFunction = staticEmFunction
   }

}

---------
-- RUN --
---------
momentApp:run()
