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

   local hasCflFrac = tbl.cflFrac and true or false
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0 -- for internal use
   self.errEps = tbl.errEps and tbl.errEps or 1.0e-8 -- residual norm
   self.maxSteps = tbl.maxSteps and tbl.maxSteps or 10000 -- max number of steps

   self.fact = tbl.factor and tbl.factor or 100 -- time-step factor over explicit
   self.extraStages = tbl.extraStages and tbl.extraStages or 1 -- extra stages
   
   -- extrapolate every these many steps
   self.extrapolateInterval = tbl.extrapolateInterval and tbl.extrapolateInterval or self.maxSteps+1 

   -- one of 'RKL1', 'richard2'
   self.stepper = tbl.stepper and tbl.stepper or 'RKL1'

   if self.stepper == "richard2" then

      local ndim = self.onGrid:ndim()
      
      self.fact = 1.0
      local L2 = 0
      for d = 1, ndim do
	 local len = self.onGrid:upper(d)-self.onGrid:lower(d)
	 L2 = L2 + 1/len^2
      end

      local L1 = math.sqrt(L2)

      -- Empirically derived bounds. I am not sure if this is the
      -- right thing to do, but for now this seems to work (AHH
      -- 11/11/2020)
      if hasCflFrac == false then
	 if ndim == 1 then
	    cflFrac = 1.5
	 else
	    cflFrac = 2^ndim
	 end
      end

      -- some adjustments to L1 to account for RDG scheme. See note
      -- above
      if ndim == 2 then L1 = L1/math.sqrt(2) end
      if ndim == 3 then L1 = 2*L1/math.sqrt(3) end      
      
      self.richardNu = 2*math.pi*L1 -- kmin
      print("richardNu", self.richardNu)
   end

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
   -- FIX-THIS: sync across procs
   
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
   -- FIX-THIS: sync across procs
   
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

-- Takes fIn and fDiff0 (which is calcRHS on fIn) and computes fOut
function IterPoisson:sts(dt, fIn, fDiff0, fOut, fact)
   local numStages = self:calcNumStages(fact, self.extraStages)
   local fDiff = self.fDiff
   local fJ, fJ1, fJ2 = self.fJ, self.fJ1, self.fJ2

   -- stage 1 (fDiff0 is already computed in main loop).
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

-- Takes fIn1, fIn and fDiff0 (which is calcRHS on fIn) and computes fOut
-- (This function uses central difference for the df/dt term)
function IterPoisson:richard2(dt, fIn1, fIn, fDiff0, fOut)
   local nu = self.richardNu

   -- compute various factors
   local fcp = (1/dt^2+nu/dt)
   local fcm = (1/dt^2-nu/dt)

   -- note that fDiff0 is already computed in main-loop. So here, all
   -- we do is compute the iteration
   fOut:combine(2/dt^2/fcp, fIn, -fcm/fcp, fIn1, 1/fcp, fDiff0)
   self:applyBc(fOut)
end

-- Takes fIn1, fIn and fDiff0 (which is calcRHS on fIn) and computes fOut
-- (This function uses a first-order fwd Euler for the df/dt term)
function IterPoisson:richard2a(dt, fIn1, fIn, fDiff0, fOut)
   local nu = self.richardNu

   -- compute various factors
   local fcp = (1/dt^2+2*nu/dt)
   local fc2 = (2/dt^2+2*nu/dt)

   -- note that fDiff0 is already computed in main-loop. So here, all
   -- we do is compute the iteration
   fOut:combine(fc2/fcp, fIn, -1/dt^2/fcp, fIn1, 1/fcp, fDiff0)
   self:applyBc(fOut)
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

   -- compute mean integrated source
   local gridVol = 1.0
   for d = 1, grid:ndim() do
      gridVol = gridVol*(grid:upper(d)-grid:lower(d))
   end
   local srcInt = self:integrateField(src)/gridVol

   -- we need to adjust sources when all directions are
   -- periodic. (FIX-THIS when not doing periodic BC)
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

   local omegaCFL = 0.0 -- compute maximum CFL frequency
   if self.stepper == "RKL1" then
      for d = 1, grid:ndim() do
	 omegaCFL = omegaCFL + 1/grid:dx(d)^2
      end
   else
      local dx1 = 0.0
      for d = 1, grid:ndim() do
	 omegaCFL = omegaCFL + 1/grid:dx(d)
      end

   end
   
   local isDone = false
   local numPrevStored = 0 -- flag to check if we are ready to extrapolate
   local errE1, errE2 = 1e10, 1e10 -- errors for use in extrapolation

   local numStages = 1
   if self.stepper == "RKL1" then
      numStages = self:calcNumStages(self.fact, self.extraStages)
   end

   if self.verbose then
      if self.stepper == "RKL1" then
	 print(string.format(" Number of stages per-step are %d", numStages))
      else
	 print(string.format(" Using Richarson second-order iteration"))
      end
   end

   local f, fNew = self.f, self.fNew
   local fE1, fE2 = self.fE1, self.fE2
   local errHist, extraHist = self.errHist, self.extraHist

   local tmStart = Time.clock()
   
   f:copy(fOut)
   local cfl = self.cfl
   local step = 0
   while not isDone do
      local dt = cfl/omegaCFL

      -- compute RHS here so we can compute the residual norm properly
      self:calcRHS(f, self.fDiff0)
      local err = self:l2norm(self.fDiff0)/srcL2
      errHist:appendData(numStages*step, { err } )      

      if self.verbose then
	 print(string.format("  Step %d, dt = %g. Res. norm = %g", step, dt, err))
      end

      if err < self.errEps or step>=self.maxSteps then
   	 isDone = true
	 break
      end      

      -- take one iteration
      if self.stepper == "RKL1" then
	 self:sts(dt, f, self.fDiff0, fNew, self.fact)
      else
	 self:richard2(dt, self.fJ1, f, self.fDiff0, fNew)
      end
      self.fJ1:copy(f) -- for richard2 scheme
      f:copy(fNew)

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
