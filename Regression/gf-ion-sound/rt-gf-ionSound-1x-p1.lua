-- Gkyl -------------------------------------------------------------------
--
-- Gyrofluid ion sound wave test case
--
---------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").Gyrofluid()

q_e     = -1.0
q_i     = 1.0
m_e     = 1.0
m_i     = 1.0

TiTe    = 1.0
kpar    = 0.5
XL, XU  = -math.pi/kpar, math.pi/kpar
Ti0     = 1.0
Te0     = Ti0/TiTe
n0      = 1.0
ni0     = n0
ne0     = n0
B0      = 1.0

vtElc = math.sqrt(Te0/m_e)
vtIon = math.sqrt(Ti0/m_i)

-- Beer & Hammett 3+1 closure model parameters.
beta_par = (32. - 9.*math.pi)/(3.*math.pi - 8.)
D_par    = 2.*math.sqrt(math.pi)/(3.*math.pi - 8.)
D_perp   = math.sqrt(math.pi)/2.

kappaParIon  = ni0*vtIon*(3.+beta_par)/(math.sqrt(3)*D_par*kpar)
kappaPerpIon = ni0*vtIon/(math.sqrt(2)*D_perp*kpar)
kappaParElc  = ne0*vtElc*(3.+beta_par)/(math.sqrt(3)*D_par*kpar)
kappaPerpElc = ne0*vtElc/(math.sqrt(2)*D_perp*kpar)

-- Geometry related parameters.
R0 = 1.0
r0 = 0.0

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = 15.0,              -- End time.
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
   cflFrac     = 0.90,

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = q_i,  mass = m_i,
      kappaPar = kappaParIon,  kappaPerp = kappaPerpIon,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = function (t, xn)
            local x = xn[1]
            local k = kpar
            local alpha = 0.01
            local perturb = alpha*math.cos(k*x)
            return ni0*(1+perturb)
         end,
         parallelTemperature = function (t, xn) return Ti0 end,
         perpendicularTemperature = function (t, xn) return Ti0 end,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","upar","Tpar","Tperp","ppar","pperp"},
   },

   -- Gyrokinetic electronss.
   elc = Plasma.Species {
      charge = q_e,  mass = m_e,
      kappaPar = kappaParElc,  kappaPerp = kappaPerpElc,
      -- Initial conditions.
      init = Plasma.GyrofluidProjection {
         density = function (t, xn) return ne0 end,
         parallelTemperature = function (t, xn) return Te0 end,
         perpendicularTemperature = function (t, xn) return Te0 end,
      },
      evolve = true, -- Evolve species?
      diagnostics = {"intMom","intM0","intM1","intM2","M2flow","upar","Tpar","Tperp","ppar","pperp"},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
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
