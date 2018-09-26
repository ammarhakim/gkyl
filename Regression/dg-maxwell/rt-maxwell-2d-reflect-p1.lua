-- Gkyl ------------------------------------------------------------------------
local Vlasov = require "App.VlasovOnCartGrid"

-- Physical parameters
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1.0

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 3.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {-1.0, -1.0}, -- configuration space lower left
   upper = {1.0, 1.0}, -- configuration space upper right
   cells = {48, 48}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- field solver
   field = Vlasov.EmField {
      epsilon0 = epsilon0, mu0 = mu0,
      elcErrorSpeedFactor = 0,
      mgnErrorSpeedFactor = 0,
      
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local rad2 = x^2 + y^2
	 return 0.0, 0.0, math.exp(-25*rad2), 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Vlasov.EmField.bcReflect, Vlasov.EmField.bcReflect },
      bcy = { Vlasov.EmField.bcReflect, Vlasov.EmField.bcReflect },
   },   
   
}
-- run application
vlasovApp:run()
