-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()

mass = 1.0
grav = 2.0
T0 = 1.0
n0 = 0.5

L = 2.5

endTime = 0.5

App = Vlasov.App {
   tEnd = 1.0,
   nFrame = 10,
   lower = {-L},
   upper = {L},
   cells = {16},
   basis = 'serendipity',
   polyOrder = 2,
   timeStepper = 'rk3',
   decompCuts = {1},
   useShared = false,
   periodicDirs = {1},   

   neut = Vlasov.Species {
      charge = 0.0, mass = 1.0,   

      lower = {-5.0*T0},
      upper = {5.0*T0},
      cells = {16},   

      init = Vlasov.MaxwellianProjection {
         density = function (t, xn)
            local x = xn[1]
            return n0*math.exp(-(x*x)/(2*0.5*0.5))
         end,
         driftSpeed = {0.0},
         temperature = T0,
         isInit = true,
      },   
   
      evolve = true,
      diagnostics = {"M0", "M1i", "M2"},   
      
      vlasovExtForceFunc = function(t, xn)
         return -5.0, 0.0, 0.0
      end,
   },   

}

App:run()
