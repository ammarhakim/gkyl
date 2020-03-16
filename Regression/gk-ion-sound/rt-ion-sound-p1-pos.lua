-- GK ion sound wave test case
-- to check result, use 
-- pgkyl -f ion-sound_phi2_ -f ion-sound_phi2-correct.h5 log plot

local Plasma = require "App.PlasmaOnCartGrid"

TiTe = 1.0
knumber = 0.5
XL, XU = -math.pi/knumber, math.pi/knumber
Ti0 = 1.0
Te0 = Ti0/TiTe
n0 = 1.0
ni0 = n0
ne0 = n0
B0 = 1.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 15, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {16}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- gyrokinetic ions
   ion = Plasma.GkSpecies {
      charge = 1.0,
      mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {32},
      -- initial conditions
      -- specify background so that we can plot perturbed distribution and moments
      initBackground = {"maxwellian",
              density = function (t, xn)
                 return ni0
              end,
              temperature = function (t, xn)
                 return Ti0
              end,
             },
      init = {"maxwellian",
              density = function (t, xn)
                 local x, v = xn[1], xn[2]
                 local k = knumber
                 local alpha = .8
                 local perturb = alpha*math.cos(k*x)
                 return ni0*(1+perturb)
              end,
              temperature = function (t, xn)
                 return Ti0
              end,
             },
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkM2", perturbed = false},
      diagnosticIntegratedMoments = {"intDelL2", "intDelPosL2"},
      positivity = true,
      positivityDiffuse = true,
   },

   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = -1.0,
      mass = 1.0,
      temp = Te0,
      -- initial conditions.. use ion background so that background is exactly neutral
      init = function (t, xn)
         return ne0
      end,
      evolve = false, -- evolve species?
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve field?
      kperp2 = 0.0,
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
