-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Left/right state for shock.
local nl, ul, pl = 1.0, 0.0, 1.0
local nr, ur, pr = 0.125, 0.0, 0.1

local vthl = math.sqrt(pl/nl)
local vthr = math.sqrt(pr/nr)

local K = 0.1

local maxwellianL = function(t, z)
   return nl/math.sqrt(2*math.pi*vthl^2)*math.exp(-(z[2]-ul)^2/(2*vthl^2))
end
local maxwellianR = function(t, z)
   return nr/math.sqrt(2*math.pi*vthr^2)*math.exp(-(z[2]-ur)^2/(2*vthr^2))
end

sim = Plasma.App {
   logToFile   = false,
   --cflFrac     = 0.1,
   tEnd        = 0.11,          -- End time.
   nFrame      = 1,             -- Number of frames to write.
   lower       = {0.0},         -- Configuration space lower left.
   upper       = {1.0},         -- Configuration space upper right.
   cells       = {128},         -- Configuration space cells.
   basis       = "serendipity", -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,             -- Polynomial order.
   timeStepper = "rk3",         -- One of "rk2", "rk3" or "rk3s4".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space
   periodicDirs = {},  -- Periodic directions.

   -- Electrons.
   neut = Plasma.Species {
      --nDiagnosticFrame = 2,
      charge = 0.0,  mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = { 6.0},
      cells = {32},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
	 density     = 0.5,--0.5*(nl+nr),
	 driftSpeed  = {0.5*(nl+nr)},
	 temperature = 0.5*(vthl^2+vthr^2),
      },
      -- reservoir = Plasma.MaxwellianProjection {
      --    density = function (t, z)
      --       if math.abs(z[1]) < 0.5 then return nl else return nr end
      --    end,
      --    driftSpeed = function (t, z)
      --       if math.abs(z[1]) < 0.5 then return {ul} else return {ur} end
      --    end,
      --    temperature = function (t, z)
      --       if math.abs(z[1]) < 0.5 then return vthl^2 else return vthr^2 end
      --    end,
      --    isReservoir = true,
      -- },
      bcx = { maxwellianL, maxwellianR },
      -- Evolve species?
      evolve = true,
      diagnostics = { "M0", "M1i", "M2" },
     
      -- Collisions.
      coll = Plasma.BGKCollisions {
         collideWith = {"neut"},
	 frequencies = {vthl/K},
      },
   },
}
-- Run application.
sim:run()
