-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell()

-- left/right state for shock
local nl, ul, pl = 1.0, 0.0, 1.0
local nr, ur, pr = 0.125, 0.0, 0.1

local vthl = math.sqrt(pl/nl)
local vthr = math.sqrt(pr/nr)

local K = 0.1


sim = Plasma.App {
   logToFile = false,
   --cflFrac = 0.1,
   tEnd = 0.11, -- end time
   nFrame = 1, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {128}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   neut = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {32},
      decompCuts = {1},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
	 density = 0.5,--0.5*(nl+nr),
	 driftSpeed = {0.5*(nl+nr)},
	 temperature = 0.5*(vthl^2+vthr^2),
      },
      reservoir = Plasma.MaxwellianProjection {
	 density = function (t, z)
	    if math.abs(z[1]) < 0.5 then return nl else return nr end
	 end,
	 driftSpeed = function (t, z)
	    if math.abs(z[1]) < 0.5 then return {ul} else return {ur} end
	 end,
	 temperature = function (t, z)
	    if math.abs(z[1]) < 0.5 then return vthl^2 else return vthr^2 end
	 end,
	 isReservoir = true,
      },
      bcx = { Plasma.Species.bcReservoir,
	      Plasma.Species.bcReservoir },
       -- evolve species?
      evolve = true,
      --evolveCollisions = false,
      -- diagnostic moments
      diagnosticMoments = { "M0", "M1i", "M2" },
     
      -- collisions
      coll = Plasma.BGKCollisions {
         collideWith = {"neut"},
	 frequencies = {vthl/K},
      },
   },
}
-- run application
sim:run()
