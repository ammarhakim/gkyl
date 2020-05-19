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

local function createField(grid, basis, vComp, onDevice, ghosts)
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
      createDeviceCopy = onDevice,
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
   local pLower = {0.0, -6.0}
   local pUpper = {1.0,  6.0}
   local pN     = {8, 16}

   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({pLower[1]}, {pUpper[1]}, {pN[1]})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distF      = createField(phaseGrid, phaseBasis, true, false)
   local numDensity = createField(confGrid, confBasis, true, false)

   -- Updater to initialize distribution function.
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         -- Using a projected Maxwellian does not give exactly 1,
         -- so would need to use assert_close instead of assert_equal. 
--         local x, v = xn[1], xn[2]
--         local n0   = 1.0
--         local u    = 0.0
--         local vt   = 0.8
--         local k    = 2.0
--         local den = n0 --*math.cos(2.0*math.pi*k*x)
--	 return (den/math.sqrt(2.0*math.pi*vt^2))
--               *math.exp( -((v-u)^2)/(2.0*(vt^2)) )

         -- Can instead project a function whose zeroth moment is exactly 1.
         return 1/(phaseGrid:upper(2)-phaseGrid:lower(2))
      end)
   project:advance(0.0, {}, {distF})

   -- Copy distribution function to device.
   local err = distF:copyHostToDevice()

   -- Moment updater.
   local calcNumDensity = distFmoment(phaseGrid,phaseBasis,confBasis,"M0",true)
   calcNumDensity:advance(0.0, {distF}, {numDensity})

   err = cudaRunTime.DeviceSynchronize()

   -- Copy moment from device to host.
   numDensity:copyDeviceToHost()

   numDensity:write("numDensity.bp",0.0)
   local momIdxr = numDensity:genIndexer()
   local nItr    = numDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking M0 moment, 1x1v")

   err = cudaRunTime.DeviceReset()
   
end

local polyOrderMax = 1
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
