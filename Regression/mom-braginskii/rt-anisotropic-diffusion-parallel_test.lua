local Grid       = require "Grid"
local Basis      = require "Basis"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"

local kappaPara = 0.01
local kappaPerp = 0.0
local timeStepper = 'symmetric-cell-center'

local lower        = { 1./2.,  1./2.}
local upper        = {-1./2., -1./2.}
local cells        = {128, 128}
local periodicDirs = {1, 2}

local cfl = 0.9
local tStart = 0
local tEnd = 3
local nFrames = 2

local tFrame = (tEnd-tStart)/nFrames

local dx = (upper[1]-lower[1]) / cells[1]
local dy = (upper[2]-lower[2]) / cells[2]
local kappa = math.max(kappaPara, kappaPerp)
local tmp = kappa / (dx)^2 + kappa / (dy)^2
local initDt = cfl * 0.5 / tmp

local function calcT (t, xn)
   local x, y = xn[1], xn[2]

   local Pi   = math.pi
   local Lx   = {upper[1] - lower[1], upper[2] - lower[2]}
   local kNum = {2*Pi/Lx[1], 2*Pi/Lx[2]}

   return math.sin(kNum[1]*x) * math.sin(kNum[2]*y)
end

local function initEmf (t, xn)
   local x, y = xn[1], xn[2]
   return 0, 0, 0, 1, 1, 0, 0, 0
end

local grid = Grid.RectCart {
   lower        = lower,
   upper        = upper,
   cells        = cells,
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim=grid:ndim(), polyOrder=0 }

local function createField(grid, basis, numComponents)
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = numComponents,
      ghost         = {2, 2},
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      },
   }
   fld:clear(0.0)
   return fld
end

local anisotropicDiffusion = Updater.AnisotropicDiffusion {
   onGrid      = grid,
   kappaPara   = kappaPara,
   kappaPerp   = kappaPerp,
   timeStepper = timeStepper,
   cfl         = cfl,
}

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function (t, xn) return 1.0 end, -- Set a dummy function for now.
   onGhosts = true,
}

-- Manually synchronize corner cells. This should really be handled by
-- CartField automatically.
local Lin = require "Lib.Linalg"
local syncCorners = function (q)
   local Nx, Ny = cells[1], cells[2]

   local idxFrom = Lin.IntVec(grid:ndim())
   local idxTo = Lin.IntVec(grid:ndim())
   local qIdxr = q:genIndexer()
   local qPtrFrom = q:get(1)
   local qPtrTo = q:get(1)

   for _,idxToFrom in ipairs({
      {   0,    0, Nx, Ny},
      {   0, Ny+1, Nx,  1},
      {Nx+1, Ny+1,  1,  1},
      {Nx+1,    0,  1, Ny},
   }) do
      idxTo[1] = idxToFrom[1]
      idxTo[2] = idxToFrom[2]
      idxFrom[1] = idxToFrom[3]
      idxFrom[2] = idxToFrom[4]
      q:fill(qIdxr(idxFrom), qPtrFrom)
      q:fill(qIdxr(idxTo), qPtrTo)
      for c=1,q:numComponents() do
         qPtrTo[c] = qPtrFrom[c]
      end
   end

   -- Fill supposedly-unused corner cells with NaNs to make sure that their
   -- values are not used in the actual update.
   local qPtr = q:get(1)
   for _,idx in ipairs({
         {  -1,   -1}, {  -1,    0}, {   0,   -1},
         {  -1, Ny+2}, {  -1, Nx+1}, {   0, Ny+2},
         {Nx+2, Ny+2}, {Nx+1, Ny+2}, {Nx+2, Nx+1},
         {Nx+2,   -1}, {Nx+1,   -1}, {Nx+2,    0},
   }) do
      idx[1] = idx[1]
      idx[2] = idx[2]
      q:fill(qIdxr(idx), qPtr)
      for c=1,q:numComponents() do
         qPtr[c] = 0/0
      end
   end
end

local temp = createField(grid, basis, 1)  -- temperature
local emf = createField(grid, basis, 8)  -- E, B, phiE, phiB
local buf = createField(grid, basis, 4)  -- grad(temp) or q, div(q)

project._isFirst = true
project:setFunc(calcT)
project:advance(tStart, {}, {temp})
temp:write("temp_0.bp", tStart, 0)
temp:sync()

project._isFirst = true  -- so that we can reuse 'project'
project:setFunc(initEmf)
project:advance(tStart, {}, {emf})
emf:sync()

local tCurr = tStart
local myDt = initDt
local frame = 1
local step = 1
while tCurr<tEnd do
   if tCurr+myDt > tEnd then myDt = tEnd-tCurr end
   print(string.format('step=%-4d, tCurr=%-10g, myDt=%-10g', step, tCurr, myDt))

   anisotropicDiffusion:setDtAndCflRate(myDt)
   local status, dtSuggested = anisotropicDiffusion:advance(
      tCurr, {buf}, {temp, emf})
   temp:sync()
   syncCorners(temp)

   if status then
      step = step + 1
      tCurr = tCurr + myDt

      if tCurr >= frame*tFrame or
         math.abs(tCurr-frame*tFrame) < 1e-10 then
         temp:write("temp_"..frame..".bp", tCurr, frame)
         frame = frame + 1
      end
   else
      print(" ** Time step "..myDt.." too large! Will retake with dt "..
            dtSuggested..".")
   end
   myDt = dtSuggested
end
