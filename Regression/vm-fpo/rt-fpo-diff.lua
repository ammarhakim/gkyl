local fpo = require "Proto.Fpo"

local kx, ky, kz = 1.0, 1.0, 1.0

sim = fpo {
   --cflFrac = 0.1,
   polyOrder = 1,
   tEnd = 1e-2,
   nFrames = 1,
   fixedDt = 5e-4,
   cells = {8, 8, 8},
   lower = {-math.pi/kx, -math.pi/ky, -math.pi/kz},
   upper = {math.pi/kx, math.pi/ky, math.pi/kz},

   periodicDirs = {1,2,3},
   init = function (t, z)
      return math.sin(kx*z[1]+ky*z[2]+kz*z[3])
   end,
   dynamicPotentials = false,
   fixedH = function (t, z)
      return 0.0
   end,
   fixedG = function (t, z)
      return z[1]^2 + z[2]^2 + z[3]^2 + z[1]*z[2] + z[1]*z[3] + z[2]*z[3]
   end,
   writeDiagnostics = false,
}

sim()
