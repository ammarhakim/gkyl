-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

knumber      = 0.5    -- Eave-number.
elVTerm      = 0.2    -- Electron thermal velocity.
vDrift       = 1.0    -- Drift velocity.
perturbation = 1.0e-6 -- Distribution function perturbation.

tFinal  = 8*4.0
nFrames = 1

xLower = {-math.pi/knumber}
xUpper = { math.pi/knumber}
Lx     = {xUpper[1]-xLower[1]}

vLower0 = {-6.0}
vUpper0 = { 6.0}
vCells0 = {32}
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

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd        = tFinal,             -- End time.
   nFrame      = nFrames,            -- Number of output frames.
   lower       = xLower, -- Configuration space lower left.
   upper       = xUpper, -- Configuration space upper right.
   cells       = {64},               -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- cuts in each configuration direction.
   useShared  = false, -- if to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Electrons.
   elc = Plasma.Species {
      charge = -1.0,  mass = 1.0,
      -- Velocity space grid.
      lower = vLower,
      upper = vUpper,
      cells = vCells,
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local alpha = perturbation
	 local k = knumber
	 local vt = elVTerm
	 
	 local fv = 1/math.sqrt(8*math.pi*vt^2)*(math.exp(-(v-vDrift )^2/(2*vt^2))+math.exp(-(v+vDrift)^2/(2*vt^2)))
	 return (1+alpha*math.cos(k*x))*fv
      end,
      vFlux  = "upwind",  -- Use upwind fluxes in velocity space.
      mask = maskOut,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i", "intM0", "M2", "intM2Flow", "intM1i", "intM2" },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0,  mu0 = 1.0,
      init = function (t, xn)
	 local alpha = perturbation
	 local k = knumber
	 return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
