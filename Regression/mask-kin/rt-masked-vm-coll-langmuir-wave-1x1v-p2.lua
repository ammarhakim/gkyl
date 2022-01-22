-- Gkyl ------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Maxwellian in 1x1v, you may not need this but it made the example file more compact
local function maxwellian1D(n, vx, ux, mass, temp)
   local v2 = (vx - ux)^2
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-mass*v2/(2*temp))
end

-- normalization parameters, shouldn't need to adjust
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

nuElc = 0.10 --sample electron collision frequency

-- Breaking up the initial condition into regions
n1 = 1.0 -- number density in region 1
 
Te1 = 1.0 -- electron temperature in region 1

-- derived parameters
vtElc1 = math.sqrt(Te1/elcMass)
-- plasma frequency and Debye length
wpe1 = math.sqrt(elcCharge^2*n1/(epsilon0*elcMass))
lambdaD1 = vtElc1/wpe1

-- parameters for perturbation
knumber = 0.5/lambdaD1
perturbation = 1e-4

vLower0 = {-6.0*vtElc1}
vUpper0 = { 6.0*vtElc1}
vCells0 = {64}
dv = {(vUpper0[1]-vLower0[1])/vCells0[1]}

extraCells = {6}
vLower = {vLower0[1]-(extraCells[1]/2)*dv[1]}
vUpper = {vUpper0[1]+(extraCells[1]/2)*dv[1]}
vCells = {vCells0[1]+extraCells[1]}

maskOut = function(t, xn)
   local x, v = xn[1], xn[2]
   local vMin = vLower0[1]
   local vMax = vUpper0[1]
   if v < vMin or v > vMax then
      return -1
   else
      return 1
   end
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 1.0/wpe1, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   cflFrac = 0.1,
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions
   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,

   -- electrons
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = vLower,
      upper = vUpper,
      cells = vCells,
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local alpha = perturbation
         local k = knumber

	 local fv = maxwellian1D(n1, v, 0.0, elcMass, Te1) 
	 return (1+alpha*math.cos(k*x))*fv
      end,
      mask = maskOut,
      evolve = true, -- evolve species?
      -- write out density, flow, total energy, and heat flux moments
      diagnostics = { "M0", "M1i", "M2", "intM0", "intM1i", "intM2" },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
