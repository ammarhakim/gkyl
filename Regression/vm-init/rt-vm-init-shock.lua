-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

gasGamma = 3.0 -- gas adiabatic constant
LX = 10.0 -- length of domain
mfp = 0.5 -- mean-free path

-- left state
M = 2.0 -- mach number for flow
nl = 1.0 -- left density
pl = 1.0 -- left pressure
ul = M*math.sqrt(gasGamma*pl/nl)

-- right state
nr = nl*(gasGamma+1)*M^2/((gasGamma-1)*M^2+2)
pr = pl*(2*gasGamma*M^2/(gasGamma+1) - (gasGamma-1)/(gasGamma+1))
ur = ul*(nl/nr)

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

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 0.5*LX/vThermal_l, -- end time
   nFrame = 10, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {LX}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   neutNoFix = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal},
      upper = {10.0*vThermal},
      cells = {8},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5*LX then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, u, vt, v)
      end,
      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" }
   },

   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal},
      upper = {10.0*vThermal},
      cells = {8},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local n = nl
	    if xn[1]>0.5*LX then
	       n = nr
	    end
	    return n
	 end,
	 driftSpeed = function (t, xn)
	    local u = ul
	    if xn[1]>0.5*LX then
	       u = ur
	    end
	    return { u }
	 end,
         temperature = function (t, xn)
	    local T = vThermal_l^2
	    if xn[1]>0.5*LX then
	       T = vThermal_r^2
	    end
	    return T
	 end,
         exactLagFixM012 = true,
      },
      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" }
   },   
}
