local math = require("sci.math").generic
local diff = require("sci.diff")

local function Bmag(t, xn)
   x, y = xn[1], xn[2]
   return 2
end

local function Bgrad(xn)
   local function BmagUnpack(...)
      local xn = {...}
      return Bmag(0, xn)
   end
   local df2 = diff.derivativef(BmagUnpack, #xn)
   local y, dx, dy, dz = df2(unpack(xn))
   return dx, dy, dz
end

local function bcurvX(t, xn)
   local bgradX, bgradY, bgradZ = Bgrad(xn)
   return -bgradY/Bmag(t,xn)^2
end
local function bcurvY(t, xn)
   local bgradX, bgradY, bgradZ = Bgrad(xn)
   return bgradX/Bmag(t,xn)^2
end

print(bcurvX(0, {4,1,1}))
print(bcurvY(0, {4,1,1}))
