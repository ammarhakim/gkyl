-- Gkyl ------------------------------------------------------------------------
--
-- Test the 2D ConstDiffusion equation object to solve a diffusion equation
--    u_t = D_x*u_{xx}+D_y*u_{yy}.
-- with constant diffusion coefficient D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Dependencies.
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Basis      = require "Basis"
local Eq         = require "Eq.ConstDiffusion"
local Updater    = require "Updater"


local pOrder       = 1                  -- Polynomial order.
local basisName    = "Ser"              -- Type of polynomial basis.
local numCells     = {8, 8}             -- Number of cells.
local lower        = {-1./2., -1./2.}   -- Lower domain boundaries.
local upper        = { 1./2.,  1./2.}   -- Upper domain boundaries.

local periodicDirs = {1,2}              -- Periodic directions.

local diffCoeff    = {0.01, 0.01}       -- Diffusion coefficient.

local tStart       = 0.0                -- Start time.
local tEnd         = 1.5                -- End time.
local nFrames      = 1                  -- Number of frames to write out.
local cflNum       = 1.0/(2*pOrder+1)   -- CFL factor. 

-- .................. end of user inputs (MAYBE) .................... --

local initDt     = 0.00976563/384 --tEnd-tStart -- initial time-step
local frame      = 1
local step       = 1
local tCurr      = tStart
local myDt       = initDt
local frameT     = tEnd/nFrames
local frameCount = 1

local function createGrid(lo, up, nCells, pDirs)
   local gridOut = Grid.RectCart {
      lower        = lo,
      upper        = up,
      cells        = nCells,
      periodicDirs = pDirs,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, numComp, vComp)
   numComp = numComp or basis:numBasis()
   vComp   = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = numComp*vComp,
      ghost         = {1, 1},
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      }
   }
   fld:clear(0.0)
   return fld
end


-- Grids and basis.
local phaseGrid  = createGrid(lower, upper, numCells, periodicDirs)
local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basisName)

-- Fields.
local distf    = createField(phaseGrid, phaseBasis)
local distf1   = createField(phaseGrid, phaseBasis)
local distfNew = createField(phaseGrid, phaseBasis)
local distfA   = createField(phaseGrid, phaseBasis)

-- Constant diffusion equation object.
local constDiffusionCalc = Eq {
   coefficient = diffCoeff,
   basis       = phaseBasis,
}

-- Updater to solve the (parabolic) diffusion equation.
local constDiffusionSlvr = Updater.HyperDisCont {
   onGrid   = phaseGrid,
   basis    = phaseBasis,
   cfl      = cflNum,
   equation = constDiffusionCalc,
}

-- Updater to project a function onto the basis.
local project = Updater.ProjectOnBasis {
   onGrid   = phaseGrid,
   basis    = phaseBasis,
   evaluate = function (t, xn) return 1.0 end -- Set the function with setFunc.
}

-- Initial condition to apply, centered at L/2 where L=1.0.
local function fInitial(xn, lower, upper)
   local x, y = xn[1], xn[2]
   local Pi   = math.pi
   local Lx   = {upper[1] - lower[1], upper[2] - lower[2]}
   local kNum = {2*Pi/Lx[1], 2*Pi/Lx[2]}
   -- Single sine mode.
   return math.sin(kNum[1]*x)*math.cos(kNum[2]*y)
end

-- Analytic answer to diffusing a sine wave.
function fAnalytic(t, xn, lower, upper, Dcoeff)
   local x, y = xn[1], xn[2]
   local Pi   = math.pi
   local Lx   = {upper[1] - lower[1], upper[2] - lower[2]}
   local kNum = {2*Pi/Lx[1], 2*Pi/Lx[2]}
   -- Single sine mode.
   return math.sin(kNum[1]*x)*math.cos(kNum[2]*y)*math.exp(-(Dcoeff[1]*kNum[1]^2+Dcoeff[2]*kNum[2]^2)*t)
end

-- Project initial condition.
project:setFunc(function(t,xn) return fInitial(xn, lower, upper) end)
project:advance(0.0, {}, {distf})
distf:write("distf_0.bp", tStart, 0)  -- Write out initial condition.
distf:sync()

-- Project the analytic solution.
project:setFunc(function(t,xn) return fAnalytic(t, xn, lower, upper, diffCoeff) end)
project:advance(tCurr, {}, {distfA})
distfA:write("distfA_0.bp", tCurr)

local cflRateByCell = createField(phaseGrid, phaseBasis, 1)
local dt = GKYL_MAX_DOUBLE

-- Function to take a single forward-euler time-step.
local function forwardEuler(tCurr, dt, fIn, fOut)
   local calcCflFlag = false
   local dtSuggested
   if dtIn == nil then calcCflFlag = true end
   cflRateByCell:clear(0.0)

   constDiffusionSlvr:setDtAndCflRate(dt, cflRateByCell)
   constDiffusionSlvr:advance(tCurr, {fIn}, {fOut})

   if calcCflFlag then
      dtSuggested = tEnd - tCurr + 1e-20

      -- Calculate suggested dt from cflRateByCell.
      -- Loop over local region.
      local grid = phaseGrid
      local cfl  = cflNum
      dt = GKYL_MAX_DOUBLE
      local localRange  = cflRateByCell:localRange()
      local cflRatePtr  = cflRateByCell:get(1)
      local cflRateIdxr = cflRateByCell:genIndexer()

      for idx in localRange:colMajorIter() do
         -- Calculate min dt from cflRateByCell.
         cflRateByCell:fill(cflRateIdxr(idx), cflRatePtr)
         dt = math.min(dt, cfl/cflRatePtr:data()[0])
      end
      dtSuggested = math.min(dt, dtSuggested)
   else
      dtSuggested = dtIn -- From argument list.
   end
   -- Take forward Euler step.
   -- NOTE: order of these arguments matters... fOut must come before fIn
   fOut:combine(dtSuggested, fOut, 1.0, fIn)
   fOut:sync()
   return dtSuggested
end

-- Function to advance solution using SSP-RK3 scheme.
local function rk3(tCurr, dt)
   local status, dtSuggested
   -- RK stage 1
   dt = forwardEuler(tCurr, dt, distf, distf1)

   -- RK stage 2
   forwardEuler(tCurr+dt, dt, distf1, distfNew)
   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)

   -- RK stage 3
   forwardEuler(tCurr+dt/2, dt, distf1, distfNew)
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   distf:copy(distf1)

   return true, dt
end

local tmSimStart = Time.clock()
-- Main simulation loop.
while true do
   -- If needed adjust dt to hit tEnd exactly.
   if tCurr+myDt > tEnd then myDt = tEnd-tCurr end
   local status, dtSuggested = rk3(tCurr, myDt)   -- Take a time-step.
   if status then   -- Check if step was successful.
      tCurr = tCurr + myDt
      -- We comment this out to take more steps and get more frames.
      myDt = dtSuggested
      step = step + 1
      if tCurr >= frameCount*frameT or math.abs(tCurr-frameCount*frameT) < 1e-10 then
         distf:write(string.format("distf_%d.bp",frameCount), tCurr, frameCount)
         project:advance(tCurr, {}, {distfA})
         distfA:write(string.format("distfA_%d.bp",frameCount), tCurr, frameCount)
         frameCount = frameCount+1
      end
      if (tCurr >= tEnd) then break end
   else
      print(string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
      myDt = dtSuggested
   end
end -- End of time-step loop.
local tmSimEnd = Time.clock()
print(string.format(" Number of steps taken: %i", step))

print(string.format(" Total time for simulation %g\n", tmSimEnd - tmSimStart))
