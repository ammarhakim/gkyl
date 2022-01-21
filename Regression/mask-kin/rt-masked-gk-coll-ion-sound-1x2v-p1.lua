-- Gkyl -------------------------------------------------------------------
--
-- GK ion acoustic wave damping.
-- Time normalized to the ion gyrofrequency, omega_ci 
-- Distances normalized to the ion sound gyroradius, rho_s=c_s/omega_ci.
-- Temperatures normalized to the electron temperature, T_e.
-- Velocities normalized to c_s=sqrt(T_e/m_i).
--
---------------------------------------------------------------------------

local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic()

tau     = 1.0            -- Ion temperature to electron temperature ratio.
Te      = 1.0            -- Electron temperature.
mIon    = 1.0            -- Ion mass
mElc    = 1.0/1836.16    -- Electron mass.
n0      = 1.0            -- Background density.
B0      = 1.0            -- Background magnetic field magnitude.
nuIon   = 2.0            -- Ion-ion collision frequency.
knumber = 0.5            -- Wave number.

-- Lower and upper bound of configuration space.
xLower, xUpper = -math.pi/knumber, math.pi/knumber

Ti      = Te*tau     -- Ion temperature.
nIon    = n0         -- Background ion density.
nElc    = n0         -- Background electron density.

vtIon   = math.sqrt(Ti)                -- Ion thermal speed. 
vtElc   = math.sqrt((mIon/mElc)*Te)    -- Electron thermal speed. 

tFinal  = .50        -- Final simulation time.
nFrames = 1

vLower  = {-6.0*vtIon, 0.0}
vUpper  = { 6.0*vtIon, mIon*((6.0*vtIon)^2)/(2.0*B0)}
vCells  = {60, 12}
print("vLower = ", vLower[1], vLower[2])
print("vUpper = ", vUpper[1], vUpper[2])
print("vCells = ", vCells[1], vCells[2])

local function round(num, numDecimalPlaces)
  local mult = 10^(numDecimalPlaces or 0)
  return math.floor(num * mult + 0.5) / mult
end

local vMinMask = {-5.0*vtIon, 0.0}
local vMaxMask = { 5.0*vtIon, mIon*((5.0*vtIon)^2)/(2.0*B0)} 
-- Find the indices for the cell whose boundaries are closes to the mask limits.
local iMin, iMax = {-1,-1}, {-1,-1}
local dv = {(vUpper[1]-vLower[1])/vCells[1],(vUpper[2]-vLower[2])/vCells[2]}
for d = 1, 2 do
   for i = 1, vCells[d] do
      local vLo, vUp = vLower[d]+(i-1)*dv[d], vLower[d]+i*dv[d]
      if iMin[d]==-1 then
         if math.abs(vLo-vMinMask[d]) < 1.e-12 then iMin[d] = i
         elseif math.abs(vUp-vMinMask[d]) < 1.e-12 then iMin[d] = i+1
         elseif vMinMask[d] > vLo and vMinMask[d] <= vUp then
            iMin[d] = i+round((vMinMask[d]-vLo)/dv[d])
         end
      end
      if iMax[d]==-1 then
         if math.abs(vLo-vMaxMask[d]) < 1.e-12 then iMax[d] = i
         elseif math.abs(vUp-vMaxMask[d]) < 1.e-12 then iMax[d] = i+1
         elseif vMaxMask[d] > vLo and vMaxMask[d] <= vUp then
            iMax[d] = i+round((vMaxMask[d]-vLo)/dv[d])
         end
      end
      if iMin[d]>-1 and iMax[d]>-1 then break end
   end
end
print("dv = ",dv[1],dv[2])
print("Mask details:")
print("  vpar limits = ", vLower[1]+(iMin[1]-1)*dv[1], vLower[1]+(iMax[1]-1)*dv[1]) 
print("  mu limits =   ", vLower[2]+(iMin[2]-1)*dv[2], vLower[2]+(iMax[2]-1)*dv[2])
print("  iMin = ",iMin[1],iMin[2])
print("  iMax = ",iMax[1],iMax[2])

ionMask = function(t, xn)
   local x, vpar, mu = xn[1], xn[2], xn[3]
   -- Mask vpar and mu space.
   if vpar < vLower[1]+(iMin[1]-1)*dv[1] or vpar > vLower[1]+(iMax[1]-1)*dv[1]
      or mu < vLower[2]+(iMin[2]-1)*dv[2] or mu > vLower[2]+(iMax[2]-1)*dv[2] then
   -- Mask vpar space only.
--   if vpar < vLower[1]+(iMin[1]-1)*dv[1] or vpar > vLower[1]+(iMax[1]-1)*dv[1] then
   -- Mask mu space only.
--   if mu < vLower[2]+(iMin[2]-1)*dv[2] or mu > vLower[2]+(iMax[2]-1)*dv[2] then
   -- No masking.
--   if vpar < vLower[1] or vpar > vUpper[1] then
      return -1
   else
      return 1
   end
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = tFinal,           -- End time.
   nFrame      = nFrames,                -- Number of output frames.
   lower       = {xLower},         -- Configuration space lower left.
   upper       = {xUpper},         -- Configuration space upper right.
   cells       = {8},              -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   cflFrac     = 1.0,

   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 1./(10.*nFrames),

   -- Decomposition for configuration space.
   decompCuts = {1},    -- Cuts in each configuration direction.
   useShared  = false,   -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},   -- periodic directions.

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = 1.0,  mass = mIon,
      -- Velocity space grid.
      lower = vLower,
      upper = vUpper,
      cells = vCells,
      -- Initial conditions.
      -- Specify background so that we can plot perturbed distribution and moments.
--      background = Plasma.MaxwellianProjection{
--         density     = function (t, xn) return nIon end,
--         temperature = function (t, xn) return Ti end,
--      },
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
            local x, v, mu = xn[1], xn[2], xn[3]
            local k        = knumber
            local alpha    = 0.01
            local perturb  = alpha*math.cos(k*x)
            return nIon*(1.0+perturb)
         end,
         temperature = function (t, xn) return Ti end,
      },
      coll = Plasma.LBOCollisions{
         collideWith = { 'ion' },
         frequencies = { nuIon },
      }, 
      evolve = true, -- Evolve species?
      mask = ionMask,
      diagnostics = {"M0", "M2","intM0","intM1","intM2"},
   },

   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = -1.0,  mass = mElc,
      temp   = Te,
      -- Initial conditions.. use ion background so that background is exactly neutral.
      init = function (t, xn) return nElc end,
      evolve = false, -- Evolve species?
   },

   -- Field solver.
   field = Plasma.Field {
      evolve  = true, -- Evolve field?
      kperpSq = 0.0,
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         return B0
      end,
      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
