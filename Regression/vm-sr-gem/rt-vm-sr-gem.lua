-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permiability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- Speed of light.

n0 = 1.0 -- Initial reference number density.

ionMass = 1.0 -- Positron mass.
ionCharge = 1.0 -- Positron charge.
elcMass = 1.0 -- Electron mass.
elcCharge = -1.0 -- Electron charge.

wpe = math.sqrt(n0*elcCharge^2/(epsilon0*elcMass))
wpi = math.sqrt(n0*ionCharge^2/(epsilon0*ionMass))
de = lightSpeed/wpe
di = lightSpeed/wpi

-- sigma = 1.0 -> magnetic field energy = rest-mass energy density
B0 = math.sqrt(2.0)

w0 = 1.0*di
psi0 = 0.1*B0*di
guide1 = 0.1*B0
guide2 = 0.1*B0

-- Inputs can be used for asymmetric reconnection, but just do symmetric for now
b1 = B0
b2 = b1
n1 = n0
n2 = n1
-- non-relativistic temperature to start
T_i1 = 0.04 
T_e1 = 0.04 
T_i2 = T_i1
T_e2 = T_e1

Lx = 20.0*math.pi*di
Ly = 20.0*math.pi*di
Nx = 16
Ny = 16
Ncx = 1
Ncy = 1

-- Maxwell-Juttner in 2x3v
local function maxwelljuttner3V(n, px, py, pz, ux, uy, uz, T, K_2)
  local gamma = 1.0/math.sqrt(1 - ux*ux - uy*uy - uz*uz);
  return n/(4*math.pi*T*K_2)*math.exp(-(gamma/T)*(math.sqrt(1 + px*px + py*py + pz*pz) - ux*px - uy*py - uz*pz));
end

-- initial conditions
function init(x,y)
   local tanh = math.tanh
   local cosh = math.cosh
   local sinh = math.sinh
   local cos = math.cos
   local sin = math.sin
   local sqrt = math.sqrt

   local me = elcMass
   local mi = ionMass
   local qe = elcCharge
   local qi = ionCharge
   
   local Pi = math.pi
   local _2pi = 2.0*Pi
   local _4pi = 2.0*_2pi

   local b1x = 0.5*(b2+b1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(b2-b1)
   local b1z = 0.5*(guide2-guide1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(guide2+guide1)
   local b1y = 0.0

   local Tvali = 0.5*(T_i2-T_i1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_i2+T_i1)
   local Tvale = 0.5*(T_e2-T_e1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_e2+T_e1)
   local n = (0.5*(b1^2 - b1x^2) +0.5*(guide1^2 - b1z^2) + n1*(T_i1+T_e1))/(Tvali+Tvale)

   local Bx = b1x - psi0*_4pi/Ly*sin(_2pi*x/Lx)*sin(_4pi*y/Ly)
   local By = b1y + psi0*_2pi/Lx*cos(_2pi*x/Lx)*(1-cos(_4pi*y/Ly))
   local Bz = b1z

   -- n = n0
   -- J = curl B

   local TeFrac = Tvale/(Tvale + Tvali)
   local TiFrac = Tvali/(Tvale + Tvali)

   local Jx = 0.5*(guide2-guide1)/w0*((1.0/cosh((y-Ly*.25)/w0))^2 - (1.0/cosh((y-Ly*.75)/w0))^2 + (1.0/cosh((y-Ly*1.25)/w0))^2 - (1.0/cosh((y+Ly*.25)/w0))^2)
   local Jy = 0.
   local Jz  = -0.5*(b2+b1)/w0*((1.0/cosh((y-Ly*.25)/w0))^2 - (1.0/cosh((y-Ly*.75)/w0))^2 + (1.0/cosh((y-Ly*1.25)/w0))^2 - (1.0/cosh((y+Ly*.25)/w0))^2) - psi0*sin(_2pi*x/Lx)*((_2pi/Lx)^2*(1 - cos(_4pi*y/Ly)) + (_4pi/Ly)^2*cos(_4pi*y/Ly))
   
   local uxe = Jx*TeFrac/(qe*n) 
   local uxi = Jx*TiFrac/(qi*n) 
   local uye = Jy*TeFrac/(qe*n) 
   local uyi = Jy*TiFrac/(qi*n) 
   local uze = Jz*TeFrac/(qe*n)  
   local uzi = Jz*TiFrac/(qi*n) 

   return n, uxe, uxi, uye, uyi, uze, uzi, Bx, By, Bz
end

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd        = 0.01,              -- End time.
   nFrame      = 1,                  -- Number of output frames.
   lower       = {0.0, 0.0}, -- Configuration space lower left.
   upper       = {Lx, Ly},  -- Configuration space upper right.
   cells       = {Nx, Ny},              -- Configuration space cells.
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {Ncx, Ncy},   -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space
   periodicDirs = {1, 2}, -- Periodic directions.

   -- Electrons.
   elc = Plasma.GenSpecies.VlasovSR {
      charge = elcCharge,  mass = elcMass,
      -- Velocity space grid.
      lower = {-2.0, -2.0, -2.0},
      upper = {2.0, 2.0, 2.0},
      cells = {8, 8, 8},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, px, py, pz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12
         local n, uxe, uxi, uye, uyi, uze, uzi, Bx, By, Bz = init(x,y)
         local fv = maxwelljuttner3V(n, px, py, pz, uxe, uye, uze, T_e1, K_2)
         return fv
      end,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i" }
   },

   -- Ions (positrons)
   ion = Plasma.GenSpecies.VlasovSR {
      charge = ionCharge,  mass = ionMass,
      -- Velocity space grid.
      lower = {-2.0, -2.0, -2.0},
      upper = {2.0, 2.0, 2.0},
      cells = {8, 8, 8},
      -- Initial conditions.
      init = function (t, xn)
         local x, y, px, py, pz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12
         local n, uxe, uxi, uye, uyi, uze, uzi, Bx, By, Bz = init(x,y)
         local fv = maxwelljuttner3V(n, px, py, pz, uxi, uyi, uzi, T_i1, K_2)
         return fv
      end,
      evolve = true, -- Evolve species?
      diagnostics = { "M0", "M1i" }
   },

   -- Field solver
   field = Plasma.GenField.Maxwell {
      epsilon0 = 1.0,  mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local n, uxe, uxi, uye, uyi, uze, uzi, Bx, By, Bz = init(x,y)
         return 0.0, 0.0, 0.0, Bx, By, Bz
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
vlasovApp:run()
