-- Gkyl ------------------------------------------------------------------------
--
-- Test for CUDA kernels that perform reductions.
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
local DataStruct  = require "DataStruct"
local Grid        = require "Grid"
local Basis       = require "Basis"
local lume        = require "Lib.lume"
local cudaRunTime = require "Cuda.RunTime"
local cudaAlloc   = require "Cuda.Alloc"
local Alloc       = require "Lib.Alloc"
local ffi         = require "ffi"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

ffi.cdef [[
  void getNumBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                              int maxThreads, int &blocks, int &threads);
  void cartFieldMax(int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads, GkDeviceProp *prop,
                    GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);
]]

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

local function createField(grid, basis, deviceCopy, ghosts, isP0, vComp)
   vComp = vComp or 1
   if ghost then
      ghostCells = {1, 1}
   else
      ghostCells = {0, 0}
   end
   if isP0 then
      numBasisElements = 1
   else
      numBasisElements = basis:numBasis()*vComp
   end
   local fld = DataStruct.Field {
      onGrid           = grid,
      numComponents    = numBasisElements,
      ghost            = ghostCells,
      createDeviceCopy = deviceCopy,
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      },
   }
   return fld
end

local function nextPowerOf2(nIn)
   return math.floor(2^(math.ceil(math.log(nIn,2))))
end

local function prevPowerOf2(nIn)
   return math.floor(2^(math.floor(math.log(nIn,2))))
end


local pOrder        = 1
local basis         = "Ser"
local phaseLower    = {0.0, -6.0}
local phaseUpper    = {1.0,  6.0}
local phaseNumCells = {3, 31}
local confLower     = {   phaseLower[1]}
local confUpper     = {   phaseUpper[1]}
local confNumCells  = {phaseNumCells[1]}

local devNum, _    = cudaRunTime.GetDevice()
local err          = cudaRunTime.SetDevice(devNum)
local devProp, err = cudaRunTime.GetDeviceProperties(devNum)

local numCellsTot   = lume.reduce(phaseNumCells, function(a,b) return a*b end)
local numBlocksMAX  = 64
local numThreadsMAX = GKYL_DEFAULT_NUM_THREADS

-- Phase-space and config-space grids.
local phaseGrid = createGrid(phaseLower, phaseUpper, phaseNumCells)
local confGrid  = createGrid(confLower, confUpper, confNumCells)
-- Basis functions.
local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
-- Fields.
local p0Field = createField(phaseGrid, phaseBasis, true, false, true)
-- Initialize field to random numbers.
math.randomseed(1000*os.time())
local fldRange = p0Field:localRange()
local fldIdxr  = p0Field:genIndexer()
for idx in fldRange:rowMajorIter() do
   local fldItr = p0Field:get(fldIdxr( idx ))
   fldItr[1]    = math.random() 
end
p0Field:copyHostToDevice()
-- Get the maximum on the CPU (for reference).
maxVal = -9.0e19
for idx in fldRange:rowMajorIter() do
   local fldItr = p0Field:get(fldIdxr( idx ))
   maxVal       = math.max(maxVal,fldItr[1])
end
print(string.format(" CPU max = %f",maxVal))
print("")


print(" Operation on the GPU:")
-- Establish the number of blocks and threads to use.
local numBlocksC   = Alloc.Int(1)
local numThreadsC  = Alloc.Int(1)
ffi.C.getNumBlocksAndThreads(devProp,numCellsTot,numBlocksMAX,numThreadsMAX,numBlocksC:data(),numThreadsC:data());
local numBlocks  = numBlocksC[1]
local numThreads = numThreadsC[1]
print(string.format("   Use %d blocks of size %d (threads)",numBlocks,numThreads))

-- Allocate device memory needed for intermediate reductions.
local d_blockRed, d_intermediateRed = cudaAlloc.Double(numBlocks), cudaAlloc.Double(numBlocks)
local d_maxVal = cudaAlloc.Double(1)

ffi.C.cartFieldMax(numCellsTot,numBlocks,numThreads,numBlocksMAX,numThreadsMAX,devProp,
                   p0Field._onDevice, d_blockRed:data(), d_intermediateRed:data(), d_maxVal:data())

-- Test that the value found is correct.
local maxVal_gpu = Alloc.Double(1)
local err        = d_maxVal:copyDeviceToHost(maxVal_gpu)

print(string.format("   GPU max = %f",maxVal_gpu[1]))
print("")
print(string.format(" Error = %e",maxVal-maxVal_gpu[1]))

cudaRunTime.Free(d_blockRed)
cudaRunTime.Free(d_intermediateRed)
cudaRunTime.Free(d_maxVal)

