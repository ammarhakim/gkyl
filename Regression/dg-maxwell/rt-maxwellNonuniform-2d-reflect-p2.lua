-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Physical parameters.
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1.0

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 3.0,            -- End time.
   nFrame = 1,            -- Number of output frames.
   lower = {-1.0, -1.0},  -- Configuration space lower left.
   upper = { 1.0,  1.0},  -- Configuration space upper right.
   coordinateMap = {
      function(z) return math.sin((2.*math.pi/4.)*z) end,
      function(z) return math.sin((2.*math.pi/4.)*z) end,
   },
   cells = {24, 24},      -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder = 2,         -- Polynomial order.
   timeStepper = "rk3",   -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space
   decompCuts = {1, 1},   -- Cuts in each configuration direction.

   -- Field solver.
   field = Vlasov.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      elcErrorSpeedFactor = 0,
      mgnErrorSpeedFactor = 0,
      
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rad2 = x^2 + y^2
         return 0.0, 0.0, math.exp(-25*rad2), 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
      bcx = { Vlasov.Field.bcReflect, Vlasov.Field.bcReflect },
      bcy = { Vlasov.Field.bcReflect, Vlasov.Field.bcReflect },
   },   
   
}
-- Run application.
vlasovApp:run()
