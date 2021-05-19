-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
-- For runs in parallel use:
-- ~/gkylsoft/openmpi/bin/mpirun -n 4 ~/gkylsoft/gkyl/bin/gkyl *lua

local Plasma = (require "App.PlasmaOnCartGrid").Gyrokinetic()
--local VmPlasma = (require "App.PlasmaOnCartGrid").VlasovMaxwell()
local Constants = require "Lib.Constants"

function sech(x)
   return 2*math.exp(x)/(math.exp(2*x)+1)
end


nfac = 1.0
-- Universal parameters.
eps0 = Constants.EPSILON0
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV
qi   =  eV
mi   = Constants.PROTON_MASS   -- Hydrogen ions.
me   = Constants.ELECTRON_MASS --*1836.0
B0  = 0.5     -- Magnetic field amplitude [T].

n0  = nfac*5e18    -- Number density [1/m^3].

Te0 = 40*eV   -- Electron temperature.
Ti0 = 40*eV   -- Ion temperature.

-- Derived physical parameters
-- set kmin*rho_s = 0.2, used for Field table
cs = math.sqrt(Te0/mi)
omega_ci = eV*B0/mi
rho_s = cs/omega_ci
kmin = 0.2/rho_s

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)

-- Connection length
Lx = 40 --[m]

-- Source parameters
Ls = 25 -- [m]
S0 = 7e22*nfac -- [m^-3 s^-1] jupyter nb recommends 
P_src = 1e6
srcFac = 0.1

-- Source profiles
sourceDensity = function (t,xn)
   local x = xn[1]
   return srcFac*S0*(math.exp(-x^2/20)+ 0.001)
   -- local flr = 0.001
   -- if math.abs(x) < Ls/2 then
   --    return S0*(math.cos(math.pi*x/Ls)+flr)
   -- else
   --  return S0*flr
   -- end
end
sourceTemperatureElc = function (t,xn)
   local x = xn[1]
   return 4.0*Te0/srcFac
end
sourceTemperatureIon = function (t,xn)
   local x = xn[1]
   return 4.0*Te0/srcFac
end
sourceDensityNeut = function (t,xn)
   local x = xn[1]
   local x0 = 1
   local S0n = 8e-21*n0^2 --nfac*1e23
   return S0n
   --  if x <= 0 then
   --    return S0n*(sech((-Lx/2-x)/x0)^2+ 0.001)
   -- else       
   --    return S0n*(sech((Lx/2-x)/x0)^2 + 0.001)
   -- end
end

-- Parameters for collisions.
nuFrac = 0.1
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))

--print(Lx/cs)

sim = Plasma.App {
   logToFile = false,

   tEnd        = 1e-7,  --1000/omega_pe,    -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = {-Lx/2},            -- Configuration space lower left.
   upper       = {Lx/2}, -- Configuration space upper right.
   cells       = {32},            -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 2,                -- Polynomial order.
   cflFrac     = 1,                -- CFL "fraction". Usually 1.0.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {1},   -- Cuts in each configuration direction.
   useShared  = false, -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {}, -- Periodic directions.
      
    -- Electrons.
   elc = Plasma.Species {
      nDistFuncFrame = 10,
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {-4.0*vte, 0},
      upper = {4.0*vte, 12*me*vte^2/(2*B0)},
      cells = {16, 8},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
	    local x, vpar = xn[1], xn[2]
	    local L0 = Lx/4
	    local c_ss = math.sqrt(5/3*Te0/mi)
	    local nPeak = 1.5*n0
	    --return nPeak*(math.exp(-x^2/40)+ 0.001)
	    local perturb = 0 
	    if math.abs(x) <= L0 then
	       return nPeak*(1+math.sqrt(1-(x/L0)^2))/2*(1+perturb)
	    else
	       return nPeak/2*(1+perturb)
	    end
	 end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 0
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Te0
         end,
      },
      source = Plasma.MaxwellianProjection {
	  density = sourceDensity,
	   temperature = sourceTemperatureElc,
      }, 
      evolve = true, -- Evolve species?
      coll = Plasma.LBOCollisions {
         collideWith = {'elc', 'ion'},
         nuFrac = 0.1,
	 --frequencies = {nuElc},
      },
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
      				     "intM2" },
      diagnosticBoundaryFluxMoments = {"GkM0", "GkUpar"},
      diagnosticIntegratedBoundaryFluxMoments = {"intM0"},
      ionization = Plasma.Ionization {
      	 collideWith  = {"neut"},
      	 electrons    = "elc",
      	 neutrals     = "neut",
      	 elemCharge   = eV, 
      	 elcMass      = me,
      	 plasma       = "H",         
      },
      bcx = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},
   },

   -- Ions
   ion = Plasma.Species {
      nDistFuncFrame = 10,
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {-4.0*vti, 0},
      upper = {4.0*vti, 12*mi*vti^2/(2*B0)},
      cells = {16, 8},
      decompCuts = {1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local x, vpar = xn[1], xn[2]
	    local L0 = Lx/4
	    local c_ss = math.sqrt(5/3*Te0/mi)
	    local nPeak = 1.5*n0
	    --return nPeak*(math.exp(-x^2/40)+ 0.001)
	    local perturb = 0 
	    if math.abs(x) <= L0 then
	       return nPeak*(1+math.sqrt(1-(x/L0)^2))/2*(1+perturb)
	    else
	       return nPeak/2*(1+perturb)
	    end
	 end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 0
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Ti0
         end,
      },
      source = Plasma.MaxwellianProjection {
	  density = sourceDensity,
	  temperature = sourceTemperatureIon,
      }, 
      evolve = true,
      coll = Plasma.LBOCollisions {
         collideWith = {'ion','elc'},
         nuFrac = 0.1,
	 --frequencies = {nuIon},
      },
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
				     "intM2" },
      diagnosticBoundaryFluxMoments = {"GkM0", "GkUpar"},
      diagnosticIntegratedBoundaryFluxMoments = {"intM0"},
      ionization = Plasma.Ionization {
      	 collideWith  = {"neut"},
      	 electrons    = "elc",
      	 neutrals     = "neut",
      	 elemCharge   = eV,
      	 elcMass      = me,
      	 plasma       = "H",
      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"neut"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      	 charge = qi,
      },
      bcx = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},
   },

   neut = Plasma.Vlasov {
      nDistFuncFrame = 10,
      charge = 0.0, mass = mi,
      -- Velocity space grid
      lower = {-4.0*vti, -4.0*vti, -4.0*vti},
      upper = {4.0*vti, 4.0*vti, 4.0*vti},
      cells = {8, 8, 8},
      decompCuts = {1},
      init = Plasma.VmMaxwellianProjection {
         density = function (t, xn)
            local x, vpar = xn[1], xn[2]
            local n_n = 4*n0
	    local x0 = 0.1
	    local flr = 1e12
	    if x <= 0 then
	       return n_n*(sech((-Lx/2-x)/x0)^2) + flr
	    else	       
	       return n_n*(sech((Lx/2-x)/x0)^2) + flr
	    end
	    -- return 0.1*n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return {0,0,0}
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 10*eV
         end,
      },
      evolve = true,
      source = Plasma.VmMaxwellianProjection {
       density = sourceDensityNeut,
       driftSpeed = function (t, xn)
          return {0,0,0}
       end,
       temperature = function (t, xn)
          return 10*eV
       end,
      },
      -- coll = Plasma.VmBGKCollisions {
      --    collideWith = {'neut'},
      -- 	 frequencies = {1e6},
      -- },
      bcx = {Plasma.Vlasov.bcRecycle, Plasma.Vlasov.bcRecycle},
      recycleTemp = 10*eV,
      recycleFrac = 0.5,
      recycleIon = "ion",
      recycleTime = .5e-3,
      diagnosticMoments = { "M0", "u", "vtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1i",
      				     "intM2Flow", "intM2Thermal" },
      diagnosticBoundaryFluxMoments = {"M0"},
      diagnosticIntegratedFluxMoments = {"M0"},
      diagnosticIntegratedBoundaryFluxMoments = {"intM0"},
      ionization = Plasma.Ionization {
      	 collideWith  = {"elc"},
      	 electrons    = "elc",
      	 neutrals     = "neut",
      	 elemCharge   = eV, 
      	 elcMass      = me,
      	 plasma       = "H",         
      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"ion"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      	 charge = 0,
      },
   },
   
   -- Field solver.
   field = Plasma.Field {
      bcLowerPhi = {{ T ="N", V = 0.0}},
      bcUpperPhi = {{ T ="N", V = 0.0}},
      evolve = true,
      kperp2 = kmin*kmin,
   },

      -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- Run application.
sim:run()
