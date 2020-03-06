local Constants = require "Constants"

------------------------
-- DEFAULT PARAMETERS --
------------------------

local function setdefault(tbl, key, val)
   if tbl[key] == nil then
      tbl[key] = val
   end
   return tbl[key]
end

local function setdefaultEarth(tbl)
   setdefault(tbl, "planetRadius", 6378.1e3)
   setdefault(tbl, "planetBx0", 0)
   setdefault(tbl, "planetBy0", 0)
   setdefault(tbl, "planetBz0", -3.12e-5)
   setdefault(tbl, "rhoIn", 5e6 * Constants.PROTON_MASS)
   setdefault(tbl, "pIn", 3e-12)
   setdefault(tbl, "vxIn", 400e3)
   setdefault(tbl, "vyIn", 0)
   setdefault(tbl, "vzIn", 0)
   setdefault(tbl, "BxIn", 0)
   setdefault(tbl, "ByIn", 0)
   setdefault(tbl, "BzIn", -5e-9)

   setdefault(tbl, "rCut", 0.5 * tbl.planetRadius)
   setdefault(tbl, "xmir", 15 * tbl.planetRadius)
   setdefault(tbl, "stretch", 2)
   setdefault(tbl, "r_ramp1", 12 * tbl.planetRadius)
   setdefault(tbl, "r_ramp2", 14 * tbl.planetRadius)

   setdefault(tbl, "rInOut", tbl.planetRadius)
end

local function setdefaultMercury(tbl)
   setdefault(tbl, "planetRadius", 2439.7e3)
   setdefault(tbl, "planetBx0", 0)
   setdefault(tbl, "planetBy0", 0)
   setdefault(tbl, "planetBz0", -195e-9)
   setdefault(tbl, "rhoIn", 56e6 * Constants.PROTON_MASS)
   setdefault(tbl, "pIn", 0.19e-9)
   setdefault(tbl, "vxIn", 400e3)
   setdefault(tbl, "vyIn", 50e3)
   setdefault(tbl, "vzIn", 0)
   setdefault(tbl, "BxIn", -15.21e-9)
   setdefault(tbl, "ByIn", 8e-9)
   setdefault(tbl, "BzIn", -8.51e-9)

   setdefault(tbl, "rCut", 0.5 * tbl.planetRadius)
   setdefault(tbl, "xmir", 5 * tbl.planetRadius)
   setdefault(tbl, "stretch", 2)
   setdefault(tbl, "r_ramp1", 2 * tbl.planetRadius)
   setdefault(tbl, "r_ramp2", 2.5 * tbl.planetRadius)

   setdefault(tbl, "rInOut", tbl.planetRadius)
end

local function setdefaultGanymede(tbl)
   setdefault(tbl, "planetRadius", 2634.1e3)
   setdefault(tbl, "planetBx0", -18e-9)
   setdefault(tbl, "planetBy0", 51.8e-9)
   setdefault(tbl, "planetBz0", -716.8e-9)
   setdefault(tbl, "rhoIn", 56e6 * Constants.PROTON_MASS)
   setdefault(tbl, "pIn", 3.8e-12)
   setdefault(tbl, "vxIn", 140e3)
   setdefault(tbl, "vyIn", 0)
   setdefault(tbl, "vzIn", 0)
   setdefault(tbl, "BxIn", 0)
   setdefault(tbl, "ByIn", -6e-9)
   setdefault(tbl, "BzIn", -77e-9)

   setdefault(tbl, "rCut", 0.5 * tbl.planetRadius)
   setdefault(tbl, "xmir", 2.5 * tbl.planetRadius)
   setdefault(tbl, "stretch", 1)
   setdefault(tbl, "r_ramp1", 2 * tbl.planetRadius)
   setdefault(tbl, "r_ramp2", 2.5 * tbl.planetRadius)

   setdefault(tbl, "rInOut", tbl.planetRadius)
end

local function setdefaultObject(tbl)
   setdefault(tbl, "gasGamma", 5./3.)

   setdefault(tbl, "mu0", Constants.MU0)
   if tbl.epsilon0 == nil then
      setdefault(tbl, "lightSpeed", Constants.SPEED_OF_LIGHT / 50)
      setdefault(tbl, "epsilon0", 1 / mu0 / (tbl.lightSpeed^2))
   end

   if type(tbl.objectName) == "string" then
      if string.lower(tbl.objectName) == "earth" then
         setdefaultEarth(tbl)
      elseif string.lower(tbl.objectName) == "mercury" then
         setdefaultMercury(tbl)
      elseif string.lower(tbl.objectName) == "ganymede" then
         setdefaultGanymede(tbl)
      else
         assert(false, "Object name \""..tbl.objectName.."\" not recognized")
      end
   end

   if tbl.planetBx0 == nil or tbl.planetBy0 == nil or tbl.planetBz0 == nil then
      -- TODO make sure theta & phi are consistent with conventions
      local B0 = tonumber(tbl.planetB0)
      local theta = tonumber(tbl.planetB0theta)
      local phi = tonumber(tbl.planetB0phi)
      tbl.planetBx0 = B0 * math.sin(theta) * math.cos(phi)
      tbl.planetBy0 = B0 * math.sin(theta) * math.sin(phi)
      tbl.planetBz0 = B0 * math.cos(theta)
   end

   local nSpecies = #tbl.moments

   local mass = tbl.mass
   local totalMass = 0
   tbl.massFractions = {}
   for s = 1, nSpecies do
      totalMass = totalMass + mass[s]
   end
   for s = 1, nSpecies do
      tbl.massFractions[s] = mass[s] / totalMass
   end
end

--------------------
-- INITIALIZATION --
--------------------

local function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz, rCut)
   local xx = x - x0
   local yy = y - y0
   local zz = z - z0
   local rr = math.sqrt(xx^2 + yy^2 + zz^2)

   if rr and (rr < rCut) then
      return 0, 0, 0
   end

   local Bx = (3*xx*Dx*xx + 3*xx*Dy*yy + 3*xx*Dz*zz - Dx*rr^2) / rr^5
   local By = (3*yy*Dx*xx + 3*yy*Dy*yy + 3*yy*Dz*zz - Dy*rr^2) / rr^5
   local Bz = (3*zz*Dx*xx + 3*zz*Dy*yy + 3*zz*Dz*zz - Dz*rr^2) / rr^5

   return Bx, By, Bz
end

local function buildInitMirdip(tbl)
   -- planet parameters
   local R = tbl.planetRadius
   -- dipole streghth; B0 is B field at equator
   local Dx = tbl.planetBx0 * R ^ 3
   local Dy = tbl.planetBy0 * R ^ 3
   local Dz = tbl.planetBz0 * R ^ 3

   -- solar wind parameters
   local rhoIn = tbl.rhoIn
   local vxIn = tbl.vxIn
   local vyIn = tbl.vyIn
   local vzIn = tbl.vzIn
   local pIn = tbl.pIn

   local BxIn = tbl.BxIn
   local ByIn = tbl.ByIn
   local BzIn = tbl.BzIn

   -- mirdip setup parameters
   local rCut = tbl.rCut
   local xmir = tbl.xmir
   local stretch = tbl.stretch
   local r_ramp1  = tbl.r_ramp1
   local r_ramp2  = tbl.r_ramp2

   -- multifluid parameters
   local massFractions = tbl.massFractions
   local pressureFractions = tbl.pressureFractions
   local moments = tbl.moments
   local nSpecies = #moments
   assert(#massFractions == #pressureFractions)

   -- constants
   local gasGamma = tbl.gasGamma

   local calcRho = tbl.calcRho
   if not calcRho then
      calcRho = function(x, y, z)
         return rhoIn
      end
   end

   local calcV = tbl.calcV
   if not calcV then
      calcV = function(x, y, z)
         local xx = x > 0 and x / stretch or x
         local r = math.sqrt(xx^2 + y^2 + z^2)
         local s = (r - r_ramp1) / (r_ramp2 - r_ramp1)
         s = math.max(s, 0)
         s = math.min(s, 1)
         local vx = vxIn * s
         local vy = vyIn * s
         local vz = vzIn * s
         return vx, vy, vz
      end
   end

   local calcP = tbl.calcP
   if not calcP then
      calcP = function(x, y, z)
         return pIn
      end
   end

   local calcB0 = tbl.calcBt
   if not calcB0 then
      calcB0 = function(x, y, z)
         return dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, rCut)
      end
   end

   local calcStaticEB = tbl.calcBt
   if not calcStaticEB then
      calcStaticEB = function(x, y, z)
         local Ex, Ey, Ez = 0, 0, 0
         local Bx, By, Bz = calcB0(x, y, z)
         local phiE, phiB = 0, 0
         return Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
      end
   end

   local calcBt = tbl.calcBt
   if not calcBt then
      calcBt = function(x, y, z)
         local Bxt, Byt, Bzt = 0, 0, 0

         if x > xmir then
            local Bxd, Byd, Bzd = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, rCut)
            Bxt, Byt, Bzt = Bxt + Bxd, Byt + Byd, Bzt + Bzd

            local Bxdm, Bydm, Bzdm = dipoleB(x, y, z, 2 * xmir, 0, 0, -Dx, Dy, Dz,
                                             rCut)
            Bxt, Byt, Bzt = Bxt + Bxdm, Byt + Bydm, Bzt + Bzdm
         end

         Bxt, Byt, Bzt = Bxt + BxIn, Byt + ByIn, Bzt + BzIn

         return Bxt, Byt, Bzt
      end
   end

   local init = tbl.init
   if not init then
      init = function(x, y, z)
         local rho = calcRho(x,y,z)
         local vx,vy,vz = calcV(x,y,z)
         local p = calcP(x,y,z)

         local values = {}
         local comp = 0

         for s = 1, nSpecies do
            local rho_s = rho * massFractions[s]
            local rhovx_s = rho_s * vx
            local rhovy_s = rho_s * vy
            local rhovz_s = rho_s * vz
            local p_s = p * pressureFractions[s]

            values[comp + 1] = rho_s
            values[comp + 2] = rhovx_s
            values[comp + 3] = rhovy_s
            values[comp + 4] = rhovz_s

            local moment = moments[s]
            if moment == 5 then
               local e_s = p_s / (gasGamma - 1) +
                           0.5 * (rhovx_s^2 + rhovy_s^2 + rhovz_s^2) / rho_s
               values[comp + 5] = e_s
            elseif moment == 10 then
               values[comp + 5] = p_s +   rhovx_s^2 / rho_s -- Pxx_s
               values[comp + 6] = rhovx_s * rhovy_s / rho_s -- Pxy_s
               values[comp + 7] = rhovx_s * rhovz_s / rho_s -- Pxz_s
               values[comp + 8] = p_s +   rhovy_s^2 / rho_s -- Pyy_s
               values[comp + 9] = rhovy_s * rhovz_s / rho_s -- Pyz_s
               values[comp + 10] = p_s +   rhovz_s^2 / rho_s -- Pzz_s
            else
               assert(false)
            end
            comp = comp + moment
         end

         local Bxt, Byt, Bzt = calcBt(x, y, z)
         local Bx0, By0, Bz0 = calcB0(x, y, 1)
         local Bx1, By1, Bz1 = Bxt - Bx0, Byt - By0, Bzt - Bz0

         local Ex = - vy * Bzt + vz * Byt
         local Ey = - vz * Bxt + vx * Bzt
         local Ez = - vx * Byt + vy * Bxt

         values[comp + 1] = Ex
         values[comp + 2] = Ey
         values[comp + 3] = Ez

         values[comp + 4] = Bx1
         values[comp + 5] = By1
         values[comp + 6] = Bz1

         local phiE = 0
         local phiB = 0
         values[comp + 7] = phiE
         values[comp + 8] = phiB

         return values
      end
   end

   local initFluid = tbl.initFluid
   if not initFluid then
      initFluid = function(x, y, z, s)
         local rho = calcRho(x,y,z)
         local vx,vy,vz = calcV(x,y,z)
         local p = calcP(x,y,z)

         local rho_s = rho * massFractions[s]
         local rhovx_s = rho_s * vx
         local rhovy_s = rho_s * vy
         local rhovz_s = rho_s * vz
         local p_s = p * pressureFractions[s]

         local moment = moments[s]
         if moment == 5 then
            local e_s = p_s / (gasGamma - 1) +
                        0.5 * (rhovx_s^2 + rhovy_s^2 + rhovz_s^2) / rho_s

            return rho_s, rhovx_s, rhovy_s, rhovz_s, e_s
         elseif moment == 10 then
            local Pxx_s = p_s +   rhovx_s^2 / rho_s -- local Pxx_s
            local Pxy_s = rhovx_s * rhovy_s / rho_s -- local Pxy_s
            local Pxz_s = rhovx_s * rhovz_s / rho_s -- local Pxz_s
            local Pyy_s = p_s +   rhovy_s^2 / rho_s -- local Pyy_s
            local Pyz_s = rhovy_s * rhovz_s / rho_s -- local Pyz_s
            local Pzz_s = p_s +   rhovz_s^2 / rho_s -- local Pzz_s

            return rho_s, rhovx_s, rhovy_s, rhovz_s, Pxx_s, Pxy_s, Pxz_s, Pyy_s,
                   Pyz_s, Pzz_s
         else
            assert(false)
         end
      end
   end

   local initField = tbl.initField
   if not initField then
      initField = function(x, y, z)
         local rho = calcRho(x,y,z)
         local vx,vy,vz = calcV(x,y,z)
         local p = calcP(x,y,z)

         local Bxt, Byt, Bzt = calcBt(x, y, z)
         local Bx0, By0, Bz0 = calcB0(x, y, 1)
         local Bx1, By1, Bz1 = Bxt - Bx0, Byt - By0, Bzt - Bz0

         local Ex = - vy * Bzt + vz * Byt
         local Ey = - vz * Bxt + vx * Bzt
         local Ez = - vx * Byt + vy * Bxt

         local phiE = 0
         local phiB = 0

         return Ex, Ey, Ez, Bx1, By1, Bz1, phiE, phiB
      end
   end

   return {
      calcRho = calcRho,
      calcV = calcV,
      calcP = calcP,
      calcBt = calcBt,
      calcB0 = calcB0,
      calcStaticEB = calcStaticEB,
      init = init,
      initFluid = initFluid,
      initField = initField,
   }
end

local calcValuesIn = function(tbl)
   -- solar wind parameters
   local rhoIn = tbl.rhoIn
   local vxIn = tbl.vxIn
   local vyIn = tbl.vyIn
   local vzIn = tbl.vzIn
   local pIn = tbl.pIn

   local BxIn = tbl.BxIn
   local ByIn = tbl.ByIn
   local BzIn = tbl.BzIn

   local gasGamma = tbl.gasGamma

   -- multifluid parameters
   local massFractions = tbl.massFractions
   local pressureFractions = tbl.pressureFractions
   local moments = tbl.moments
   local nSpecies = #moments

   local valuesIn = {}

   for s = 1, nSpecies do
      local rho_s = rhoIn * massFractions[s]
      local rhovx_s = rho_s * vxIn
      local rhovy_s = rho_s * vyIn
      local rhovz_s = rho_s * vzIn
      local p_s = pIn * pressureFractions[s]

      local moment = moments[s]
      if moment == 5 then
         local e_s = p_s / (gasGamma - 1) +
                     0.5 * (rhovx_s^2 + rhovy_s^2 + rhovz_s^2) / rho_s
         valuesIn[s] = {rho_s, rhovx_s, rhovy_s, rhovz_s, e_s}
      elseif moment == 10 then
         local Pxx_s = p_s +   rhovx_s^2 / rho_s -- local Pxx_s
         local Pxy_s = rhovx_s * rhovy_s / rho_s -- local Pxy_s
         local Pxz_s = rhovx_s * rhovz_s / rho_s -- local Pxz_s
         local Pyy_s = p_s +   rhovy_s^2 / rho_s -- local Pyy_s
         local Pyz_s = rhovy_s * rhovz_s / rho_s -- local Pyz_s
         local Pzz_s = p_s +   rhovz_s^2 / rho_s -- local Pzz_s

         valuesIn[s] = {rho_s, rhovx_s, rhovy_s, rhovz_s, Pxx_s, Pxy_s, Pxz_s,
                        Pyy_s, Pyz_s, Pzz_s}

      else
         assert(false)
      end
   end

   local ExIn = -vyIn * BzIn + vzIn * ByIn
   local EyIn = -vzIn * BxIn + vxIn * BzIn
   local EzIn = -vxIn * ByIn + vyIn * BxIn
   local phiEIn = 0
   local phiBIn = 0

   valuesIn[nSpecies + 1] = {ExIn, EyIn, EzIn, BxIn, ByIn, BzIn, phiEIn, phiBIn}

   return valuesIn
end

local buildInOutFunc = function(tbl)
   return function(t, xn)
      local x, y, z = xn[1], xn[2], xn[3]
      if (x^2 + y^2 + z^2 < tbl.rInOut^2) then
         return -1
      else
         return 1
      end

   end
end

----------
-- Grid --
----------

local function calcGrid(calcDx, xlo, xup, nx)
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
      xx_nc[i + 1] = xx_nc[i] + calcDx((j + .5) / nx)
   end
   for j = 0, -sw1 + 1, -1 do
      local i = j + sw1 + 1
      xx_nc[i - 1] = xx_nc[i] - calcDx((j - .5) / nx)
   end
   local fac = 1 / xx_nc[nx + sw1 + 1]
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

local function linearRefine(xx, xx_nc, x_nc)
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
         x = x_nc[i - 1] + (x_nc[i] - x_nc[i - 1]) * (xx - xx_nc[i - 1]) / (xx_nc[i] - xx_nc[i - 1])
         break
      end
   end
   return x
end

local function calcDxFlatTop1d(x, x0, width, flatWidth, ratioDx)
   if (math.abs(x - x0) <= flatWidth) then
      return 1
   elseif (x - x0 > flatWidth) then
      local x1 = x0 + flatWidth
      return (1 / ratioDx - (1 / ratioDx - 1) * math.exp(-.5 * ((x - x1) / width) ^ 2))
   elseif (x - x0 < -flatWidth) then
      local x1 = x0 - flatWidth
      return (1 / ratioDx - (1 / ratioDx - 1) * math.exp(-.5 * ((x - x1) / width) ^ 2))
   end
end

local function buildGridFlatTop1d(xlo, xup, nx, x0, width, flatWidth, ratioDx)
   local calcDx = function(x)
      return calcDxFlatTop1d(x, (x0 - xlo) / (xup - xlo), width / (xup - xlo),
                             flatWidth / (xup - xlo), ratioDx)
   end

   local xx_nc, x_nc = calcGrid(calcDx, xlo, xup, nx)

   local coordinateMap = function(xx)
      return linearRefine(xx, xx_nc, x_nc)
   end

   return coordinateMap
end

return {
   setdefaultObject = setdefaultObject,
   buildInitMirdip = buildInitMirdip,
   buildInOutFunc = buildInOutFunc,
   calcValuesIn = calcValuesIn,
   buildGridFlatTop1d = buildGridFlatTop1d,
}
