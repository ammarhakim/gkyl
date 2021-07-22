-- Gkyl -------------------------------------------------------------------
--
-- Test the gyrofluid pitch-angle scattering operator by creating a plasma
-- with a pressure anisotropy and relaxing it.
--
---------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").Gyrofluid()

q_i     = 1.0
m_i     = 1.0

kpar    = 0.5
XL, XU  = -math.pi/kpar, math.pi/kpar
Ti0     = 1.0
n0      = 1.0
B0      = 1.0

pRat = 8./3.  -- Ratio of pPerp to pPar.

TiPar0  = 3.*Ti0/(2.*pRat+1.)
TiPerp0 = pRat*TiPar0 

vtIon = math.sqrt(Ti0/m_i)

nuPAS_ii = 0.1  -- Pitch-angle scattering rate.

-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

kappaParIon  = n0*(vtIon^2)*(3.+beta_par)/(math.sqrt(3)*D_par*vtIon*kpar)
kappaPerpIon = n0*(vtIon^2)/(math.sqrt(2)*D_perp*vtIon*kpar)

-- Geometry related parameters.
R0 = 1.0
r0 = 0.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = 60.0,              -- End time.
   nFrame = 1,                 -- Number of output frames.
   lower  = {-math.pi/kpar},   -- Configuration space lower left.
   upper  = { math.pi/kpar},   -- Configuration space upper right.
   cells  = {16},              -- Configuration space cells.
   mapc2p = function(xc)  -- MF (2021/04/20): this is only intended to give a Jacobian=1 for now.
      -- Field-aligned coordinates (x,y).
      local x, y = xc[1], xc[2]
      -- Cylindrical coordinates (R,phi).
      local phi = x/(R0+r0)
      -- Cartesian coordinates (X,Y).
      local X = R0*math.cos(phi)
      local Y = R0*math.sin(phi)
      return X, Y
   end,
   basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                  -- Polynomial order.
   timeStepper = "rk3",              -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = q_i,  mass = m_i,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = function (t, xn) return n0 end,
         parallelTemperature = function (t, xn)
            local x       = xn[1]
            local k       = kpar
            local alpha   = 0.0
            local perturb = alpha*math.cos(k*x)
            return TiPar0*(1+perturb)
         end,
         perpendicularTemperature = function (t, xn)
            local x       = xn[1]
            local k       = kpar
            local alpha   = 0.0
            local perturb = alpha*math.cos(k*x)
            return TiPerp0*(1+perturb)
         end,
      },
      closure = Plasma.HeatFlux{
         kappaPar = kappaParIon,  kappaPerp = kappaPerpIon,
      },
      coll = Plasma.PASCollisions {
         collideWith = {'ion'},
         frequencies = {nuPAS_ii},
      },
      evolve = true, -- Evolve species?
      evolveCollisionless = false,
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","Upar","Tpar","Tperp","Ppar","Pperp"},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = false, -- Evolve field?
      externalPhi = function (t, xn) return 0.0 end,
      kperpSq = 0.3,
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn) return B0 end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
