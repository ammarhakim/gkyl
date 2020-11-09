-- Gkyl ------------------------------------------------------------------------
--
-- Iterative Poisson solver
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local Eq = require "Eq.ConstDiffusion"
local HyperDisCont = require "Updater.HyperDisCont"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local xsys = require "xsys"

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)

local function mu(j) return (2*j-1)/j end
local function nu(j) return (1-j)/j end
local function mubar(s,j) return (2*j-1)/j*2/(s^2+s) end


-- Iterative Poisson solver updater object
local IterPoisson = Proto(UpdaterBase)

function IterPoisson:calcNumStages(dhdp, extraStages)
   return math.ceil(1/2*(math.sqrt(1+8*dhdp)-1)) + extraStages
end

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
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0 -- for internal use
   self.errEps = tbl.errEps and tbl.errEps or 1.0e-8 -- residual norm
   self.maxSteps = tbl.maxSteps and tbl.maxSteps or 10000 -- max number of steps

   self.fact = tbl.factor and tbl.factor or 100 -- time-step factor over explicit
   self.extraStages = tbl.extraStages and tbl.extraStages or 1 -- extra stages
   
   -- extrapolate every these many steps
   self.extrapolateInterval = tbl.extrapolateInterval and tbl.extrapolateInterval or self.maxSteps+1 

   -- one of 'RKL1', 'Richardson2'
   self.stepper = tbl.stepper and tbl.stepper or 'RKL1'

   -- flag to print internal iteration steps
   self.verbose = xsys.pickBool(tbl.verbose, false)
   
   self.cfl = self.fact*cflFrac/(2*polyOrder+1)/self.onGrid:ndim()

   local coeff = {}
   for d = 1, self.onGrid:ndim() do
      coeff[d] = 1.0
   end

   -- create updater to use in inner loop
   local constDiff = Eq {
      coefficient = coeff,
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

   -- source copy (should get rid of this)
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
   local vol = self.onGrid:cellVolume()
   local dfact = 1/math.sqrt(2^self.onGrid:ndim())
   
   local intf = 0.0
   for idxs in localRange:colMajorIter() do
      local fldItr = fld:get(indexer(idxs))
      intf = intf + fldItr[1]
   end
   -- SYNC across procs
   
   return intf*vol*dfact
end

-- compute the L2 norm
function IterPoisson:l2norm(fld)
   local localRange = fld:localRange()
   local indexer = fld:genIndexer()
   local vol = self.onGrid:cellVolume()
   local dfact = 1/2^self.onGrid:ndim()
   
   local l2 = 0.0
   for idxs in localRange:colMajorIter() do
      local fldItr = fld:get(indexer(idxs))
      for k = 1, fld:numComponents() do
	 l2 = l2 + fldItr[k]^2
      end
   end
   -- SYNC across procs
   
   return math.sqrt(l2*vol*dfact)
end

function IterPoisson:l2diff(f1, f2)
   local localRange = f1:localRange()
   local indexer = f1:genIndexer()
   local vol = self.onGrid:cellVolume()
   local dfact = 1/2^self.onGrid:ndim()

   local l2 = 0.0
   for idxs in localRange:colMajorIter() do
      local f1Itr = f1:get(indexer(idxs))
      local f2Itr = f2:get(indexer(idxs))
      
      for k = 1, f1:numComponents() do
	 l2 = l2 + (f1Itr[k]-f2Itr[k])^2
      end
   end
   return math.sqrt(l2*vol*dfact)
end

function IterPoisson:calcRHS(fIn, fOut)
   local dt, tCurr = 0.1, 1.0 -- these are ignored
   self.constDiffusionSlvr:setDtAndCflRate(dt, self.cflRateByCell)
   self.constDiffusionSlvr:advance(tCurr, {fIn}, {fOut})
   fOut:accumulate(1.0, self.src)
end

function IterPoisson:applyBc(fld)
   fld:sync()
end

function IterPoisson:sts(dt, fIn, fOut, fact)
   local numStages = self:calcNumStages(fact, self.extraStages)
   local fDiff0, fDiff = self.fDiff0, self.fDiff
   local fJ, fJ1, fJ2 = self.fJ, self.fJ1, self.fJ2

   -- we need this in each stage
   self:calcRHS(fIn, fDiff0)

   -- stage 1
   fJ2:copy(fIn)
   fJ1:combine(1.0, fIn, mubar(numStages,1)*dt, fDiff0)
   self:applyBc(fJ1)
   
   -- rest of stages
   for j = 2, numStages do
      self:calcRHS(fJ1, fDiff)
      fJ:combine(mu(j), fJ1, nu(j), fJ2, mubar(numStages,j)*dt, fDiff)
      self:applyBc(fJ)
      -- reset fields for next stage
      fJ2:copy(fJ1); fJ1:copy(fJ)
   end
   fOut:copy(fJ)
end

-- advance method
function IterPoisson:_advance(tCurr, inFld, outFld)

   local grid, basis = self.onGrid, self.basis
   local polyOrder = basis:numBasis()

   local srcIn = assert(inFld[1], "IterPoisson.advance: Must-specify an input RHS (source)")
   local fOut = assert(outFld[1], "IterPoisson.advance: Must-specify an output field")

   -- clear diagnostic data
   self.errHist:clear(); self.extraHist:clear()

   local src = self.src
   src:copy(srcIn)
   local srcInt = self:integrateField(src)/grid:cellVolume() -- mean integrated source

   -- we need to adjust sources when all directions are periodic. (FIX
   -- this when not doing periodic BC)
   do
      local localRange = src:localRange()
      local indexer = src:genIndexer()
      local dfact = math.sqrt(2^grid:ndim())
      
      for idxs in localRange:colMajorIter() do
	 local sr = src:get(indexer(idxs))
	 sr[1] = sr[1] - dfact*srcInt
      end
   end

   local srcL2 = self:l2norm(src) -- L2 norm of source

   local step = 1
   local omegaCFL = 0.0 -- compute maximum CFL frequency
   for d = 1, grid:ndim() do
      omegaCFL = omegaCFL + 1/grid:dx(d)^2
   end
   
   local isDone = false
   local numPrevStored = 0 -- flag to check if we are ready to extrapolate
   local errE1, errE2 = 1e10, 1e10 -- errors for use in extrapolation
   
   numStages = self:calcNumStages(self.fact, self.extraStages)

   if self.verbose then
      print(string.format(" Number of stages per-step are %d", numStages))
   end

   local f, fNew = self.f, self.fNew
   local fE1, fE2 = self.fE1, self.fE2
   local errHist, extraHist = self.errHist, self.extraHist

   local tmStart = Time.clock()
   
   f:copy(fOut)
   local cfl = self.cfl
   while not isDone do
      local dt = cfl/omegaCFL
      self:sts(dt, f, fNew, self.fact)
      
      local err = self:l2diff(f, fNew)
      local resNorm = err/dt/srcL2

      if self.verbose then
	 print(string.format("  Step %d, dt = %g. Error = %g (Res. norm = %g)", step, dt, err, resNorm))
      end

      f:copy(fNew)      

      if err < self.errEps or step>=self.maxSteps then
   	 isDone = true
      end

      -- check if we should store the solution for use in
      -- extrapolation
      if step % self.extrapolateInterval == 0 then
   	 fE1:copy(fE2); fE2:copy(f)
   	 errE1 = errE2; errE2 = err
   	 numPrevStored = numPrevStored + 1
   	 if numPrevStored > 1 then -- need two values to extrapolate
   	    local eps = errE2/errE1 -- extrapolation factor
   	    extraHist:appendData(numPrevStored-1, { eps } )
   	    f:combine(1.0, fE2, eps, fE2, -eps, fE1)
   	 end
      end
      errHist:appendData(numStages*step, { err } )
      step = step+1
   end

   if self.verbose then
      print(
	 string.format(
	    " IterPoisson took %g sec, %d stages", Time.clock()-tmStart, (step-1)*numStages
	 )
      )      
   end

   if step>self.maxSteps then
      assert(true, "WARNING: IterPoisson solver has not converged! Increase 'maxSteps'")
   end

   fOut:copy(f)
   
end

function IterPoisson:writeDiagnostics()
   self.errHist:write("errHist.bp", 1.0)
   self.extraHist:write("extraHist.bp", 1.0)
end

return IterPoisson
