-- Gkyl ------------------------------------------------------------------------
--
-- Iterative Poisson solver
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local DataStruct = require "DataStruct"
local Eq = require "Eq.ConstDiffusion"
local HyperDisCont = require "Updater.HyperDisCont"

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)
local function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
local function a(j) return 1-b(j) end
local function w1(s) return 4/(s^2+s-2) end
local function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end

-- For RKL2 scheme
local function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
local function nu(j) return -(j-1)/j*b(j)/b(j-2) end
local function gbar(s,j) return -a(j-1)*mubar(s,j) end

local function calcNumStagesRKL2(dhdp, extraStages) 
   return math.ceil(math.sqrt(4*dhdp+9/4) - 1/2) + extraStages
end

-- For RKL1 scheme
local function muRKL1(j) return (2*j-1)/j end
local function nuRKL1(j) return (1-j)/j end
local function mubarRKL1(s,j) return (2*j-1)/j*2/(s^2+s) end

local function calcNumStagesRKL1(dhdp, extraStages) 
   return math.ceil(1/2*(math.sqrt(1+8*dhdp)-1)) + extraStages
end

-- Iterative Poisson solver updater object
local IterPoisson = Proto(UpdaterBase)

-- constructor
function IterPoisson:init(tbl)
   IterPoisson.super.init(self, tbl) -- setup base object
   
   -- read data from input file
   self.onGrid = assert(
      tbl.onGrid, "Updater.IterPoisson: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis, "Updater.IterPoisson: Must specify basis functions to use using 'basis'")

   local polyOrder = self.basis:polyOrder()

   -- eventually determine default values with some in-built
   -- heuristics
   self.cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0 -- for internal use
   self.errEps = tbl.errEps and tbl.errEps or 1.0e-8 -- residual norm
   self.maxSteps = tbl.maxSteps and tbl.maxSteps or 10000 -- max number of steps
   self.fact = tbl.factor and tbl.factor or 100 -- time-step factor over explicit
   self.extraStages = tbl.extraStages and tbl.extraStages or 1 -- extra stages
   -- extrapolate every these many steps
   self.extrapolateInterval = tbl.extrapolateInterval and tbl.extrapolateInterval or maxSteps+1 

   -- one of 'RKL2' or 'RKL1', 'Richardson2'
   self.stepper = tbl.stepper and tbl.stepper or 'RKL1'
   
   self.cfl = self.cflFrac/(2*polyOrder+1)/self.onGrid:ndim()

   -- create updater to use in inner loop
   local constDiff = Eq {
      coefficient = {1.0, 1.0},
      basis = self.basis,
   }
   -- updater to solve (parabolic) diffusion equation: this is used in
   -- inner loop to iterate to steady-state
   self.constDiffusionSlvr = HyperDisCont {
      onGrid = self.onGrid,
      basis = self.basis,
      cfl = self.cfl,
      equation = constDiff,
   }

   -- function to allocate fields
   local function getField(numComponents)
      return DataStruct.Field {
	 onGrid = self.onGrid,
	 numComponents = numComponents,
	 ghost = {1, 1},
      }
   end

   local numBasis = self.basis:numBasis()

   -- allocate memory: need to optimize to reduce memory footprint
   self.f = getField(numBasis)
   self.fNew = getField(numBasis)
   self.fDiff0 = getField(numBasis)
   self.fDiff = getField(numBasis)
   self.fJ = getField(numBasis)
   self.fJ1 = getField(numBasis)
   self.fJ2 = getField(numBasis)

   -- for extrapolation
   self.fE1 = getField(numBasis)
   self.fE2 = getField(numBasis)

   -- source copy
   self.src = getField(numBasis)

   -- for time-stepping
   self.cflRateByCell = getField(1)

   -- error history (diagnostics)
   self.errHist = DataStruct.DynVector { numComponents = 1 }
   -- extrapolation factors (diagnostics)
   self.extraHist = DataStruct.DynVector { numComponents = 1 }
end

-- compute integral of field
function IterPoisson:integrateField(fld)
   local localRange = fld:localRange()
   local indexer = fld:genIndexer()
   local vol = grid:cellVolume()
   
   local intf = 0.0
   for idxs in localRange:colMajorIter() do
      local fldItr = fld:get(indexer(idxs))
      intf = intf + fldItr[1]/2 -- FIX
   end
   -- SYNC across procs
   
   return intf*dx*dy
end

-- compute the L2 norm
function IterPoisson:l2norm(fld)
   local localRange = fld:localRange()
   local indexer = fld:genIndexer()
   local vol = grid:cellVolume()
   
   local l2 = 0.0
   for idxs in localRange:colMajorIter() do
      local fldItr = fld:get(indexer(idxs))
      for k = 1, fld:numComponents() do
	 l2 = l2 + fldItr[k]^2/4.0 -- FIX
      end
   end
   -- SYNC across procs
   
   return math.sqrt(l2*dx*dy)
end

-- advance method
function IterPoisson:_advance(tCurr, inFld, outFld)

   -- clear diagnostic data
   self.errHist:clear(); self.extraHist:clear()
   
end

function IterPoisson:writeDiagnostics()
   self.errHist:write("errHist.bp", 1.0)
   self.extraHist:write("extraHist.bp", 1.0)
end

return IterPoisson
