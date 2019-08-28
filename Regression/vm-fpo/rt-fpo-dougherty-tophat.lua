local fpo = require "Proto.Fpo"

sim = fpo {
   cflFrac = 0.1,
   polyOrder = 2,
   tEnd = 10,
   nFrames = 1,
   cells = {12, 12},
   lower = {-6, -6},
   upper = {6, 6},
   init = function (t, z)
      if math.abs(z[1]-1) <= 1 and math.abs(z[2]-1) <= 1 then
	 return 1.0
      else
	 return 0.0
      end
   end,
   updatePotentials = false,
   doughertyPotentials = true,
   writeDiagnostics = true,
}

sim()
