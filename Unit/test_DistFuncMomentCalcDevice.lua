-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute moments on a device.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Unit        = require "Unit"
local Grid        = require "Grid"
local DataStruct  = require "DataStruct"
local Basis       = require "Basis"
local Updater     = require "Updater"
local Lin         = require "Lib.Linalg"
local cudaRunTime = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local function createGrid(lo,up,nCells)
   local gridOut = Grid.RectCart {
      lower = lo,
      upper = up,
      cells = nCells,
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

local function createField(grid, basis, deviceCopy, ghosts, vComp)
   vComp = vComp or 1
   if ghost then
      ghostCells = {1, 1}
   else
      ghostCells = {0, 0}
   end
   local fld = DataStruct.Field {
      onGrid           = grid,
      numComponents    = basis:numBasis()*vComp,
      ghost            = ghostCells,
      createDeviceCopy = deviceCopy,
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      },
   }
   return fld
end

local function createProject(grid,basis)
   local projUp = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1 end
   }
   return projUp
end

local function distFmoment(grid,pBasis,cBasis,mom,doOnDevice)
   local momUp = Updater.DistFuncMomentCalc {
      onGrid     = grid,
      phaseBasis = pBasis,
      confBasis  = cBasis,
      moment     = mom,
      onDevice   = doOnDevice
   }
   return momUp
end

local function test_moments(phaseLower, phaseUpper, phaseNumCells, 
                            confLower, confUpper, confNumCells, basis, pOrder)

   local cDim = #confLower
   local vDim = #phaseLower-cDim

   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(phaseLower, phaseUpper, phaseNumCells)
   local confGrid  = createGrid(confLower, confUpper, confNumCells)
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distF             = createField(phaseGrid, phaseBasis, true, false)
   local numDensity        = createField(confGrid, confBasis, true, false)
   local momDensity        = createField(confGrid, confBasis, true, false, vDim)
   local energyDensity     = createField(confGrid, confBasis, true, false)
   local numDensityHost    = createField(confGrid, confBasis, false, false)
   local momDensityHost    = createField(confGrid, confBasis, false, false, vDim)
   local energyDensityHost = createField(confGrid, confBasis, false, false)

   -- Updater to initialize distribution function.
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         -- Using a projected Maxwellian does not give exactly 1,
         -- so would need to use assert_close instead of assert_equal. 
         local x, v = xn[1], xn[2]
         local n0   = 1.0
         local u    = 1.0
         local vt   = 0.8
         local k    = 2.0
         local den = n0 --*math.cos(2.0*math.pi*k*x)
	 return (den/math.sqrt(2.0*math.pi*vt^2))
               *math.exp( -((v-u)^2)/(2.0*(vt^2)) )

--         -- Can instead project a function whose zeroth moment is exactly 1.
--         return 1/(phaseGrid:upper(2)-phaseGrid:lower(2))
      end)
   project:advance(0.0, {}, {distF})

   -- Copy distribution function to device.
   local err = distF:copyHostToDevice()

   local confRange    = numDensity:localRange()
   local confIdxr     = numDensity:genIndexer()

   -- ~~~~~~~~~~~~~~~~~~~ Check the number density ~~~~~~~~~~~~~~~~~~~~ --
   local calcNumDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M0",true)
   calcNumDensity:_advanceOnDevice(0.0, {distF}, {numDensity})

   err = cudaRunTime.DeviceSynchronize()

   numDensity:copyDeviceToHost()

   -- Compare against CPU calculation.
   local calcNumDensityHost = distFmoment(phaseGrid,phaseBasis,confBasis,"M0",false)
   calcNumDensityHost:_advance(0.0, {distF}, {numDensityHost})

   local m0Components = numDensity:numComponents()
   for idx in confRange:rowMajorIter() do
      local m0Itr     = numDensity:get(confIdxr( idx ))
      local m0ItrHost = numDensityHost:get(confIdxr( idx ))
      for k = 1, m0Components do
         assert_close(m0ItrHost[k], m0Itr[k], 1.e-11, string.format("Checking M0 moment, %s%dx%dvP%d",basis,cDim,vDim,pOrder))
      end
   end

--   numDensity:write("numDensity.bp",0.0)
--   numDensityHost:write("numDensityHost.bp",0.0)

   -- ~~~~~~~~~~~~~~~~~~~ Check the momentum density ~~~~~~~~~~~~~~~~~~~~ --
   local calcMomDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M1i",true)
   calcMomDensity:_advanceOnDevice(0.0, {distF}, {momDensity})

   err = cudaRunTime.DeviceSynchronize()

   momDensity:copyDeviceToHost()

   -- Compare against CPU calculation.
   local calcMomDensityHost = distFmoment(phaseGrid,phaseBasis,confBasis,"M1i",false)
   calcMomDensityHost:_advance(0.0, {distF}, {momDensityHost})

   local m1iComponents = momDensity:numComponents()
   for idx in confRange:rowMajorIter() do
      local m1iItr     = momDensity:get(confIdxr( idx ))
      local m1iItrHost = momDensityHost:get(confIdxr( idx ))
      for k = 1, m1iComponents do
         assert_close(m1iItrHost[k], m1iItr[k], 1.e-11, string.format("Checking M1i moment, %s%dx%dvP%d",basis,cDim,vDim,pOrder))
      end
   end

--   momDensity:write("momDensity.bp",0.0)
--   momDensityHost:write("momDensityHost.bp",0.0)

   -- ~~~~~~~~~~~~~~~~~~~ Check the energy density ~~~~~~~~~~~~~~~~~~~~ --
   local calcEnergyDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M2",true)
   calcEnergyDensity:_advanceOnDevice(0.0, {distF}, {energyDensity})

   err = cudaRunTime.DeviceSynchronize()

   energyDensity:copyDeviceToHost()

   -- Compare against CPU calculation.
   local calcEnergyDensityHost = distFmoment(phaseGrid,phaseBasis,confBasis,"M2",false)
   calcEnergyDensityHost:_advance(0.0, {distF}, {energyDensityHost})

   local m2Components = energyDensity:numComponents()
   for idx in confRange:rowMajorIter() do
      local m2Itr     = energyDensity:get(confIdxr( idx ))
      local m2ItrHost = energyDensityHost:get(confIdxr( idx ))
      for k = 1, m2Components do
         assert_close(m2ItrHost[k], m2Itr[k], 1.e-11, string.format("Checking M2 moment, %s%dx%dvP%d",basis,cDim,vDim,pOrder))
      end
   end

--   energyDensity:write("energyDensity.bp",0.0)
--   energyDensityHost:write("energyDensityHost.bp",0.0)
end

function test_1x1v(basis, pOrder)
   local phaseLower    = {0.0, -6.0}
   local phaseUpper    = {1.0,  6.0}
   local phaseNumCells = {8, 32}
   local confLower     = {   phaseLower[1]}
   local confUpper     = {   phaseUpper[1]}
   local confNumCells  = {phaseNumCells[1]}
   test_moments(phaseLower, phaseUpper, phaseNumCells, confLower, confUpper, confNumCells, basis, pOrder)
end

function test_1x2v(basis, pOrder)
   local phaseLower    = {0.0, -6.0, -6.0}
   local phaseUpper    = {1.0,  6.0,  6.0}
   local phaseNumCells = {8, 32, 32}
   local confLower     = {   phaseLower[1]}
   local confUpper     = {   phaseUpper[1]}
   local confNumCells  = {phaseNumCells[1]}
   test_moments(phaseLower, phaseUpper, phaseNumCells, confLower, confUpper, confNumCells, basis, pOrder)
end

function test_1x3v(basis, pOrder)
   local phaseLower    = {0.0, -6.0, -6.0, -6.0}
   local phaseUpper    = {1.0,  6.0,  6.0,  6.0}
   local phaseNumCells = {8, 32, 32, 32}
   local confLower     = {   phaseLower[1]}
   local confUpper     = {   phaseUpper[1]}
   local confNumCells  = {phaseNumCells[1]}
   test_moments(phaseLower, phaseUpper, phaseNumCells, confLower, confUpper, confNumCells, basis, pOrder)
end

function test_2x2v(basis, pOrder)
   local phaseLower    = {0.0, 0.0, -6.0, -6.0}
   local phaseUpper    = {1.0, 1.0,  6.0,  6.0}
   local phaseNumCells = {8, 8, 32, 32}
   local confLower     = {   phaseLower[1],    phaseLower[2]}
   local confUpper     = {   phaseUpper[1],    phaseUpper[2]}
   local confNumCells  = {phaseNumCells[1], phaseNumCells[2]}
   test_moments(phaseLower, phaseUpper, phaseNumCells, confLower, confUpper, confNumCells, basis, pOrder)
end

function test_2x3v(basis, pOrder)
   local phaseLower    = {0.0, 0.0, -6.0, -6.0, -6.0}
   local phaseUpper    = {1.0, 1.0,  6.0,  6.0,  6.0}
   local phaseNumCells = {8, 8, 32, 32, 32}
   local confLower     = {   phaseLower[1],    phaseLower[2]}
   local confUpper     = {   phaseUpper[1],    phaseUpper[2]}
   local confNumCells  = {phaseNumCells[1], phaseNumCells[2]}
   test_moments(phaseLower, phaseUpper, phaseNumCells, confLower, confUpper, confNumCells, basis, pOrder)
end

local polyOrderMax = 2
local basisType    = "Ser"

for polyOrder = 1, polyOrderMax do
   test_1x1v( basisType, polyOrder )
   test_1x2v( basisType, polyOrder )
   test_1x3v( basisType, polyOrder )
   test_2x2v( basisType, polyOrder )
   test_2x3v( basisType, polyOrder )
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
