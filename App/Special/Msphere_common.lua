

local function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz)
   local xx = (x - x0)
   local yy = (y - y0)
   local zz = (z - z0)
   local rr = math.sqrt(xx^2+yy^2+zz^2)
   local Bx = (3*xx*Dx*xx + 3*xx*Dy*yy + 3*xx*Dz*zz - Dx*rr^2) / rr^5
   local By = (3*yy*Dx*xx + 3*yy*Dy*yy + 3*yy*Dz*zz - Dy*rr^2) / rr^5
   local Bz = (3*zz*Dx*xx + 3*zz*Dy*yy + 3*zz*Dz*zz - Dz*rr^2) / rr^5
   return Bx, By, Bz
end

local function mirdip_dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz, cut, xmir)
   if cut and x < xmir then
      return 0, 0, 0
   end
   return dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz)
end

function mirdip_rho(x,y,z)
   local rho = rho_in
   return rho
end

function mirdip_pressure(x,y,z)
   local p = p_in
   return p
end

local function mirdip_v(x, y, z)
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

local function mirdip_staticB(x, y, z, xmir, r1)
   local r2 = x^2+y^2+z^2
   if (r2 < r1^2) then
      return 0, 0, 0
   end
   local Bxs, Bys, Bzs = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = mirdip_dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, false, xmir)
   Bxs, Bys, Bzs = Bxs+Bxd0, Bys+Byd0, Bzs+Bzd0
   -- Bxs, Bys, Bzs = Bxs+Bx_in, Bys+By_in, Bzs+Bz_in
   return Bxs, Bys, Bzs
end

local function mirdip_totalB(x, y, z, xmir, r1)
   local r2 = x^2+y^2+z^2
   if (r2 < r1^2) then
      return 0, 0, 0
   end
   local Bxt, Byt, Bzt = 0, 0, 0
   local Bxd0, Byd0, Bzd0 = mirdip_dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz, true, xmir)
   Bxt, Byt, Bzt = Bxt+Bxd0, Byt+Byd0, Bzt+Bzd0
   local Bxdm, Bydm, Bzdm = mirdip_dipoleB(x, y, z, 2*xmir, 0, 0, -Dx, Dy, Dz, true, xmir)
   Bxt, Byt, Bzt = Bxt+Bxdm, Byt+Bydm, Bzt+Bzdm
   Bxt, Byt, Bzt = Bxt+Bx_in, Byt+By_in, Bzt+Bz_in
   return Bxt, Byt, Bzt
end

local function mirdip_init5m(x, y, z, xmir, r1)
   local rho = mirdip_rho(x,y,z)
   local vx,vy,vz = mirdip_v(x,y,z)
   local p = mirdip_pressure(x,y,z)

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

   local Bxt, Byt, Bzt = mirdip_totalB(x, y, z, xmir, r1)
   local Bxs, Bys, Bzs = mirdip_staticB(x, y, z, xmir, r1)
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
   local Bxs, Bys, Bzs = mirdip_staticB(x, y, z, xmir, r1)
   return Exs, Eys, Ezs, Bxs, Bys, Bzs, 0, 0
end

function setInOutField(x, y, z)
   if (x^2 + y^2 + z^2 < r0^2) then
      return -1
   else
      return 1
   end
end

return {
   dipoleB = dipoleB,
   mirdip_staticB = mirdip_staticB,
   mirdip_init5m = mirdip_init5m,
}
