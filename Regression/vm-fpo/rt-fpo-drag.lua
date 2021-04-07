local fpo = require "Proto.Fpo"

local n = 1.0
local ux, uy, uz = 0.0, 0.0, 0.0
local vth = 1.0

local function maxwellian(vx, vy, vz, n, ux, uy, uz, vth) 
   return n/math.sqrt(2*math.pi*vth^2)^3
      * math.exp(-(vx-ux)^2/(2*vth^2))
      * math.exp(-(vy-uy)^2/(2*vth^2))
      * math.exp(-(vz-uz)^2/(2*vth^2))
end

sim = fpo {
   --cflFrac = 0.1,
   polyOrder = 1,
   tEnd = 12.0,
   nFrames = 1,
   fixedDt = 1e-1,
   cells = {8, 8, 8},
   lower = {-6.0, -6.0, -6.0},
   upper = {6.0, 6.0, 6.0},
   periodicDirs = {1,2,3},
   init = function (t, z)
      return maxwellian(z[1], z[2], z[3], n, ux, uy, uz, vth)
      --return math.sin(2*math.pi*z[1]/6)
   end,
   dynamicPotentials = false,
   fixedH = function (t, z)
      return z[1] + z[2] + z[3]
   end,
   fixedG = function (t, z)
      return 0--z[1]^2 + z[1]*z[2] + z[2]^2 + z[3]^2
   end,
   writeDiagnostics = false,
}

sim()
