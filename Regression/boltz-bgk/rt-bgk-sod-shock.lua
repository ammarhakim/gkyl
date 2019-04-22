-- Gkyl --------------------------------------------------------------
-- BGK Relexation test -----------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- left/right state for shock
nl, ul, pl = 1.0, 0.0, 1.0
nr, ur, pr = 0.125, 0.0, 0.1

Lx = 1.0 -- domain size
mfp = Lx/100 -- mean-free path

-- thermal velocity to give same energy as in fluid internal energy
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- use left state as reference
nu = vThermal/mfp -- collision frequency

VL, VU = -6.0*vThermal, 6.0*vThermal

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {Lx}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"
   cflFrac = 0.9,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- neutrals
   neut = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal},
      upper = {6.0*vThermal},
      cells = {16},
      decompCuts = {1},

      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5 then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, u, vt, v)
      end,

      evolveCollisionless = true,
      evolveCollisions = true,
      -- collisions
      bgk = Plasma.BgkCollisions {
	 collideWith = {"neut"},
	 frequencies = {nu},
      },

      bcx = { Plasma.Species.bcOpen, Plasma.Species.bcOpen },

      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
}
-- run application
sim:run()
