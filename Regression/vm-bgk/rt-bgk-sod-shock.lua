-- Gkyl --------------------------------------------------------------
--
-- BGK relaxation test
--
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Left/right state for shock.
nl, ul, pl = 1.0, 0.0, 1.0
nr, ur, pr = 0.125, 0.0, 0.1

Lx  = 1.0    -- Domain size.
mfp = Lx/100 -- Mean-free path.

-- Thermal velocity to give same energy as in fluid internal energy.
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l     -- Use left state as reference.
nu       = vThermal/mfp   -- collision frequency.

VL, VU = -6.0*vThermal, 6.0*vThermal

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd        = 0.1,             -- End time.
   nFrame      = 1,               -- Number of frames to write.
   lower       = {0.0},           -- Configuration space lower left.
   upper       = {Lx},            -- Configuration space upper right.
   cells       = {64},            -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 0.9,

   -- Decomposition for configuration space.
   decompCuts = {1},     -- Cuts in each configuration direction.
   useShared  = false,   -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {}, -- Periodic directions.

   -- Neutrals
   neut = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0*vThermal},
      upper = { 6.0*vThermal},
      cells = {16},

      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5 then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, u, vt, v)
      end,

      evolveCollisionless = true,
      evolveCollisions    = true,
      -- Collisions.
      bgk = Plasma.BGKCollisions {
	 collideWith = {"neut"},
	 frequencies = {nu},
      },

      bcx = { Plasma.OpenBC{}, Plasma.OpenBC{} },

      diagnostics = { "M0", "M1i", "M2", "M3i", "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
}
-- run application
sim:run()
