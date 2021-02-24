-- Gkyl ------------------------------------------------------------------------
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"
local Basis = require "Basis"
local EqVlasov = require "Eq.Vlasov"
local EqMaxwell = require "Eq.PerfMaxwell"
local Updater = require "Updater"
local ffi = require "ffi"

-- Constants
chargeElc = -1.0
massElc = 1.0
epsilon0 = 1.0

-- R is ratio of vth to ud
ud = 0.3
R  = 0.333333333333333

-- Initial conditions
nElc1 = 0.5
nElc2 = 0.5
uxElc1 = 0.0
uyElc1 = ud
uxElc2 = 0.0
uyElc2 = -ud
TElc1 = massElc*(R*ud)^2
TElc2 = massElc*(R*ud)^2

k0 = 1.0
theta = 45.0/180.0*math.pi --0deg is pure Weibel, 90deg is pure Two Stream
kx = k0*math.cos(theta)
ky = k0*math.sin(theta)
perturb_n = 1e-8
alpha = 1.18281106421231 --ratio of E_y/E_x 

vthElc1 = math.sqrt(TElc1/massElc)
vthElc2 = math.sqrt(TElc2/massElc)

Nx = 48 -- Number of configuration space point in x.
Ny = 48 -- Number of configuration space point in y.
Nvx = 48 -- Number of velocity space point in vx.
Nvy = 48 -- Number of velocity space point in vy.
Lx = 2*math.pi/kx -- domain size in x.
Ly = 2*math.pi/ky -- domain size in y.
Vmax = 0.9 -- Upper bound of velocity grid (3x initial drift velocity, 0.9 speed of light).

polyOrder = 2 -- polynomial order.
tEnd = 10.0 -- end time.
cflNum = 1.0/(2*polyOrder+1) -- CFL multiplicative factor.

-- NOTE: This input file assumes a fixed-size input time-step
local myDt = 0.01
local step   = 1
local tCurr  = 0.0
-- Estimate the number of time-steps from tEnd and myDt
local maxStep = tEnd/myDt

print(string.format("tEnd = %g", tEnd))
print(string.format("dt = %g", myDt))
print(string.format("Estimated number of time-steps = %g", maxStep))

-- Cartesian rectangular grid.
local phaseGrid = Grid.RectCart {
   lower = {0.0, 0.0, -Vmax, -Vmax},
   upper = {Lx, Ly, Vmax, Vmax},
   cells = {Nx, Ny, Nvx, Nvy},
   periodicDirs = {1, 2},
}
local confGrid = Grid.RectCart {
   lower = { phaseGrid:lower(1), phaseGrid:lower(2) },
   upper = { phaseGrid:upper(1), phaseGrid:upper(2) },
   cells = { phaseGrid:numCells(1), phaseGrid:numCells(2) },
   periodicDirs = {1, 2},
}
local cdim = confGrid:ndim()
local vdim = phaseGrid:ndim()-confGrid:ndim()
--  ................... Basis functions ..................... --
local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = polyOrder }
local confBasis  = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = polyOrder }

--  ......................... Fields ........................ --
-- NOTE: These first two CartFields (distribution function and em fields) are written out.
-- So, we also give them the metadata of polynomial order and basis type for easy plotting.
-- Distribution function.
local distf = DataStruct.Field {
   onGrid        = phaseGrid,
   numComponents = phaseBasis:numBasis(),
   ghost         = {1, 1},
   metaData = {
      polyOrder = phaseBasis:polyOrder(),
      basisType = phaseBasis:id()
   },
}

-- EM fields.
local em = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = 8*confBasis:numBasis(),
   ghost         = {1, 1},
   metaData = {
      polyOrder = phaseBasis:polyOrder(),
      basisType = phaseBasis:id()
   },
}

-- Momentum (1st moment).
local momDens = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = confBasis:numBasis()*vdim,
   ghost         = {1, 1},
   metaData = {
      polyOrder = phaseBasis:polyOrder(),
      basisType = phaseBasis:id()
   },
}

-- Current for coupling Maxwell's equations with Vlasov equation.
local current = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = confBasis:numBasis()*vdim,
   ghost         = {1, 1},
}
-- Total EM field (to obtain charge/mass * electromagnetic fields)
local totalEmField = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = 8*confBasis:numBasis(),
   ghost         = {1, 1},
}

-- Extra fields for RK stepping.
local distf1 = DataStruct.Field {
   onGrid        = phaseGrid,
   numComponents = phaseBasis:numBasis(),
   ghost         = {1, 1},
}
local distfNew = DataStruct.Field {
   onGrid        = phaseGrid,
   numComponents = phaseBasis:numBasis(),
   ghost         = {1, 1},
}

local em1 = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = 8*confBasis:numBasis(),
   ghost         = {1, 1},
}
local emNew = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = 8*confBasis:numBasis(),
   ghost         = {1, 1},
}

-- CFL arrays
local cflRateByCell = DataStruct.Field {
   onGrid        = phaseGrid,
   numComponents = 1,
   ghost         = {1, 1},
}
local dtGlobal = ffi.new("double[2]")

local cflRateByCellMaxwell = DataStruct.Field {
   onGrid        = confGrid,
   numComponents = 1,
   ghost         = {1, 1},
}
local dtGlobalMaxwell = ffi.new("double[2]")

-- ............................ Updaters ............................ --
VlasovCalc = EqVlasov {
      onGrid = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis = confBasis,
      charge = chargeElc,
      mass = massElc,
      hasElectricField = true,
      hasMagneticField = true,
}
local zfd = { }
for d = 1, vdim do zfd[d] = cdim+d end

MaxwellCalc = EqMaxwell {
      lightSpeed = 1.0,
      elcErrorSpeedFactor = 0.0,
      mgnErrorSpeedFactor = 0.0,
      basis = confBasis,
}

VlasovSlvr = Updater.HyperDisCont {
   onGrid             = phaseGrid,
   basis              = phaseBasis,
   cfl                = cflNum,
   equation           = VlasovCalc,
   zeroFluxDirections = zfd,
}

MaxwellSlvr = Updater.HyperDisCont {
   onGrid             = confGrid,
   basis              = confBasis,
   cfl                = cflNum,
   equation           = MaxwellCalc,
}

momDensityCalc = Updater.DistFuncMomentCalc {
   onGrid     = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis  = confBasis,
   moment     = "M1i",
}

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

-- ............................ Initial conditions ............................ -- 
function fInitial(x,y,vx,vy)
   local fv = maxwellian2D(nElc1, vx, vy, uxElc1, uyElc1, vthElc1)+maxwellian2D(nElc2, vx, vy, uxElc2, uyElc2, vthElc2)
   return (1.0+perturb_n*math.cos(kx*x+ky*y))*fv
end

function emInitial(x,y)
   local E_x = -perturb_n*math.sin(kx*x+ky*y)/(kx+ky*alpha)
   local E_y = alpha*E_x
   local B_z = kx*E_y-ky*E_x
   return E_x, E_y, 0.0, 0.0, 0.0, B_z, 0.0, 0.0
end

-- updater to initialize distribution function
local initF = Updater.ProjectOnBasis {
   onGrid   = phaseGrid,
   basis    = phaseBasis,
   onGhosts = true,
   evaluate = function (t, xn)
      return fInitial(xn[1],xn[2],xn[3],xn[4])
   end
}
local initFields = Updater.ProjectOnBasis {
   onGrid   = confGrid,
   basis    = confBasis,
   onGhosts = true,
   evaluate = function (t, xn)
      return emInitial(xn[1],xn[2])
   end
}

initF:advance(0.0, {}, {distf})
-- Compute initial momentum.
momDensityCalc:advance(0.0, {distf}, {momDens})

initFields:advance(0.0, {}, {em})

-- Output initial distribution function, momentum, and em fields.
distf:write("distf.bp", 0.0)
momDens:write("momDens.bp", 0.0)
em:write("em.bp", 0.0)

-- Copy initialized fields onto the device
distf:copyHostToDevice()
distf:devicePeriodicCopy()

em:copyHostToDevice()
em:devicePeriodicCopy()

momDens:copyHostToDevice()

-- Function to combine and accumulate forward Euler time-step.
local function combine(outIdx, a, aIdx, ...)
   local args = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   outIdx:deviceCombine(a, aIdx)
   for i = 1, nFlds do -- Accumulate rest of the fields.
      outIdx:deviceAccumulate(args[2*i-1], args[2*i])
   end
end

-- function to take a single forward-euler time-step
local function forwardEuler(tCurr, dt, fIn, emIn, fOut, emOut)
   -- Clear CFL before taking time-step.
   cflRateByCell:deviceClear(0.0)
   cflRateByCellMaxwell:deviceClear(0.0)
   -- Clear totalEmField and accumulate charge/mass electromagnetic field.
   totalEmField:deviceClear(0.0)
   local qbym = chargeElc/massElc
   totalEmField:deviceAccumulate(qbym, emIn)

   -- compute momentum from distribution function.
   momDensityCalc:advance(0.0, {fIn}, {momDens})

   -- Compute RHS of Vlasov equation.
   VlasovSlvr:setDtAndCflRate(dtGlobal[0], cflRateByCell)
   VlasovSlvr:advance(tCurr, {fIn, totalEmField}, {fOut})

   -- Compute RHS of Maxwell's equations
   MaxwellSlvr:setDtAndCflRate(dtGlobalMaxwell[0], cflRateByCellMaxwell)
   MaxwellSlvr:advance(tCurr, {emIn}, {emOut})

   -- Accumulate current onto RHS of electric field.
   -- Third input is the component offset.
   -- The offset is zero for accumulating the current onto the electric field.
   emOut:deviceAccumulateOffset(-chargeElc/epsilon0, momDens, 0)

   -- increment solution f^{n+1} = f^n + dtSuggested*fRHS
   combine(fOut, dt, fOut, 1.0, fIn)
   combine(emOut, dt, emOut, 1.0, emIn)
   -- synchronize ghost cells for periodic boundary conditions
   fOut:devicePeriodicCopy()
   emOut:devicePeriodicCopy()
end

-- function to advance solution using SSP-RK3 scheme
local function rk3(tCurr, dt)
      -- RK stage 1
      forwardEuler(tCurr, dt, distf, em, distf1, em1)

      -- RK stage 2
      forwardEuler(tCurr+dt, dt, distf1, em1, distfNew, emNew)
      distf1:deviceCombine(3.0/4.0, distf, 1.0/4.0, distfNew)
      em1:deviceCombine(3.0/4.0, em, 1.0/4.0, emNew)

      -- RK stage 3
      forwardEuler(tCurr+dt/2, dt, distf1, em1, distfNew, emNew)
      distf1:deviceCombine(1.0/3.0, distf, 2.0/3.0, distfNew)
      distf:deviceCopy(distf1)
      em1:deviceCombine(1.0/3.0, em, 2.0/3.0, emNew)
      em:deviceCopy(em1)
end

local tmSimStart = Time.clock()
-- main simulation loop
while tCurr<tEnd do
      -- take a time-step
      rk3(tCurr, myDt)
      tCurr = tCurr + myDt
      step = step + 1
      if (tCurr >= tEnd) then
         break
      end
end -- end of time-step loop
local tmSimEnd = Time.clock()

print(string.format(" Total time for simulation %g\n", tmSimEnd - tmSimStart))

-- Output the final distribution function, em fields, and momentum.
distf:copyDeviceToHost()
em:copyDeviceToHost()
momDens:copyDeviceToHost()

distf:write("distf_final.bp", 0.0)
em:write("em_final.bp", 0.0)
momDens:write("momDens_final.bp", 0.0)


