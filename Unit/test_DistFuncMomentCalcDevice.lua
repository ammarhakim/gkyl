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

function test_1x1v(pOrder, basis)
   local lower    = {0.0, -6.0}
   local upper    = {1.0,  6.0}
   local numCells = {8, 32}

   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(lower, upper, numCells)
   local confGrid  = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distF          = createField(phaseGrid, phaseBasis, true, false)
   local numDensity     = createField(confGrid, confBasis, true, false)
   local momDensity     = createField(confGrid, confBasis, true, false)
   local numDensityHost = createField(confGrid, confBasis, false, false)
   local momDensityHost = createField(confGrid, confBasis, false, false)

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

   -- ~~~~~~~~~~~~~~~~~~~ Check the number density ~~~~~~~~~~~~~~~~~~~~ --
   local calcNumDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M0",true)
   calcNumDensity:_advanceOnDevice(0.0, {distF}, {numDensity})

   err = cudaRunTime.DeviceSynchronize()

   numDensity:copyDeviceToHost()

   -- Compare against CPU calculation.
   local calcNumDensityHost = distFmoment(phaseGrid,phaseBasis,confBasis,"M0",false)
   calcNumDensityHost:_advance(0.0, {distF}, {numDensityHost})

   local numIdxr     = numDensity:genIndexer()
   local numIdxrHost = numDensityHost:genIndexer()
   for i = 1, numCells[1] do
      local numItr     = numDensity:get(numIdxr( {i} ))
      local numItrHost = numDensityHost:get(numIdxr( {i} ))
      assert_close(numItrHost[1], numItr[1], 1.e-14, "Checking M0 moment, 1x1v")
   end

   -- ~~~~~~~~~~~~~~~~~~~ Check the momber density ~~~~~~~~~~~~~~~~~~~~ --
   local calcMomDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M1i",true)
   calcMomDensity:_advanceOnDevice(0.0, {distF}, {momDensity})

   err = cudaRunTime.DeviceSynchronize()

   momDensity:copyDeviceToHost()

   -- Compare against CPU calculation.
   local calcMomDensityHost = distFmoment(phaseGrid,phaseBasis,confBasis,"M1i",false)
   calcMomDensityHost:_advance(0.0, {distF}, {momDensityHost})

   local momIdxr     = momDensity:genIndexer()
   local momIdxrHost = momDensityHost:genIndexer()
   for i = 1, numCells[1] do
      local momItr     = momDensity:get(momIdxr( {i} ))
      local momItrHost = momDensityHost:get(momIdxr( {i} ))
      assert_close(momItrHost[1], momItr[1], 1.e-14, "Checking M0 moment, 1x1v")
   end

end

local polyOrderMax = 2
local basisType    = "Ser"

for polyOrder = 1, polyOrderMax do
   test_1x1v( polyOrder, basisType )
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
