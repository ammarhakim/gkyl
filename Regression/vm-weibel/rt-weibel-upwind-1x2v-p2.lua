-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Constants
chargeElc, chargeIon = -1.0, 1.0
massElc, massIon     = 1.0, 1836.0

-- Initial conditions
nElc10, nElc20   = 0.5, 0.5
uxElc10, uyElc10 = 0.0, 0.3
uzElc10, uxElc20 = 0.0, 0.0
uyElc20, uzElc20 = -0.3, 0.0
TElc10, TElc20   = 0.01, 0.01

nIon0 = 1.0
uxIon0, uyIon0, uzIon0 = 0.0, 0.0, 0.0
TIon0 = 0.01

k0 = 0.4
perturb = 1e-3

vthElc10, vthElc20 = math.sqrt(TElc10/massElc), math.sqrt(TElc20/massElc)
vthIon0 = math.sqrt(TIon0/massIon)

knumber      = 0.5    -- Wave-number.
elVTerm      = 0.2    -- Electron thermal velocity.
vDrift       = 1.0    -- Drift velocity.
perturbation = 1.0e-6 -- Distribution function perturbation.

-- Maxwellian in 1x2v.
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd        = 50.0,             -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = { 0.0 },          -- Configuration space lower left.
   upper       = { 2*math.pi/k0 }, -- Configuration space upper right.
   cells       = {24},             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1}, -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Electrons.
   elc = Vlasov.Species {
      charge = -1.0,  mass = 1.0,
      -- Velocity space grid.
      lower = {-1.0, -1.0},
      upper = { 1.0,  1.0},
      cells = {12, 12},
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
	    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)
      end,
      evolve = true, -- Evolve species?
      vFlux  = "upwind",  -- Use upwind fluxes in velocity space.
      diagnostics = { "M0", "M1i", "M2" }
   },

   -- Field solver.
   field = Vlasov.Field {
      epsilon0 = 1.0,  mu0 = 1.0,
      init = function (t, xn)
	 local x = xn[1]
	 local Bz = perturb*math.sin(k0*x)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, Bz
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
