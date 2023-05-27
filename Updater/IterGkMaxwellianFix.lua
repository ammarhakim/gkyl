-- Gkyl ------------------------------------------------------------------------
--
-- Iterative fix for Maxwell projection
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Logger = require "Lib.Logger"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local MaxwellianOnBasis = require "Updater.MaxwellianOnBasis"

local xsys = require "xsys"

-- Iterative Poisson solver updater object
local IterGkMaxwellianFix = Proto(UpdaterBase)

-- constructor
function IterGkMaxwellianFix:init(tbl)
   IterGkMaxwellianFix.super.init(self, tbl) -- setup base object

   self.phaseGrid =
      assert(tbl.onGrid, "Updater.IterGkMaxwellianFix: Must provide phase space grid object 'onGrid'")
   self.confGrid =
      assert(tbl.confGrid, "Updater.IterGkMaxwellianFix: Must provide configuration space grid object 'confGrid'")
   self.confBasis =
      assert(tbl.confBasis, "Updater.IterGkMaxwellianFix: Must provide configuration space basis object 'confBasis'")
   self.phaseBasis =
      assert(tbl.phaseBasis, "Updater.IterGkMaxwellianFix: Must provide phase space basis object 'phaseBasis'")

   self.mass = tbl.mass
   self.bmag = tbl.bmag

   self.maxIter = tbl.maxIter
   self.relEps = tbl.relEps
   
   local phaseGrid, confGrid, confBasis, phaseBasis =
      self.phaseGrid, self.confGrid, self.confBasis, self.phaseBasis

   local mass, bmag = self.mass, self.bmag

   local polyOrder = self.phaseBasis:polyOrder()

   -- function to allocate field on conf space
   local function getField(grid, basis, numComponents)
      local f = DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis()*numComponents,
	 ghost = {1, 1},
      }
      return f
   end
   
   self.vDim = phaseGrid:ndim() - confGrid:ndim() 
   if self.vDim > 1 then
      self.vDim = 3
   end

   -- Allocate memory
   self.m0scl = getField(confGrid, confBasis, 1) 
   self.dm0 = getField(confGrid, confBasis, 1) 
   self.ddm0 = getField(confGrid, confBasis, 1) 
   self.m0_new = getField(confGrid, confBasis, 1)
   self.dm1 = getField(confGrid, confBasis, 1) 
   self.ddm1 = getField(confGrid, confBasis, 1) 
   self.m1_new = getField(confGrid, confBasis, 1)
   self.dm2 = getField(confGrid, confBasis, 1)
   self.ddm2 = getField(confGrid, confBasis, 1)
   self.m2_new = getField(confGrid, confBasis, 1)
   self.m2flow = getField(confGrid, confBasis, 1)
   self.uDrift = getField(confGrid, confBasis, 1)
   self.vtSq = getField(confGrid, confBasis, 1)
   self.vtSqTar = getField(confGrid, confBasis, 1)

   -- create updaters
   self.maxwellian = MaxwellianOnBasis {
      onGrid = phaseGrid,
      confGrid = confGrid,
      phaseBasis = phaseBasis,
      confBasis = confBasis,
      mass = mass,
   }

   -- Moment updaters.
   self.calcNumDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "GkM0",
      gkfacs = {mass, bmag},
   }
   self.calcMomDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "GkM1",
      gkfacs = {mass, bmag},
   }
   self.calcKinEnergyDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "GkM2",
      gkfacs = {mass, bmag},
   }

   -- Binary operation updaters.
   self.weakDivide = CartFieldBinOp {
      onGrid = confGrid,
      weakBasis = confBasis,
      operation = "Divide",
   }
   self.weakMultiply = CartFieldBinOp {
      onGrid = confGrid,
      weakBasis = confBasis,
      operation = "Multiply",
   }
   self.weakMultiplyConfPhase = CartFieldBinOp {
      onGrid = phaseGrid,
      weakBasis = phaseBasis,
      operation = "Multiply",
      fieldBasis = confBasis,
   }
   
   -- Add a DynVector
   self.flag = DataStruct.DynVector {
      numComponents = 3,
}
end

-- compute the L2 norm
function IterGkMaxwellianFix:l2norm(fld)
   local localRange = fld:localRange()
   local indexer = fld:genIndexer()
   local vol = self.confGrid:cellVolume()
   local dfact = 1/2^self.confGrid:ndim()

   local l2 = 0.0
   for idxs in localRange:colMajorIter() do
      local fldItr = fld:get(indexer(idxs))
      for k = 1, fld:numComponents() do
         l2 = l2 + fldItr[k]^2
      end
   end

   return math.sqrt(l2*vol*dfact)
end

function IterGkMaxwellianFix:l2diff(f1, f2)
   local localRange = f1:localRange()
   local indexer = f1:genIndexer()
   local vol = self.confGrid:cellVolume()
   local dfact = 1/2^self.confGrid:ndim()

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

function IterGkMaxwellianFix:_advance(tCurr, inFld, outFld)

   -- fM is the uncorrected Maxwellian projection
   local fM = assert(inFld[1], "IterGkMaxwellianFix.advance: Must-specify an input distribution function")
   local m0 = assert(inFld[2], "IterGkMaxwellianFix.advance: Must-specify an input number density")
   local m1 = assert(inFld[3], "IterGkMaxwellianFix.advance: Must-specify an input moment density")
   local m2 = assert(inFld[4], "IterGkMaxwellianFix.advance: Must-specify an input energy density")
   local bmag = assert(inFld[5], "IterGkMaxwellianFix.advance: Must-specify magnetic field amplitude")
      
   local fOut = assert(outFld[1], "IterGkMaxwellianFix.advance: Must-specify an output distribution function")

   local m0scl = self.m0scl
   local dm0 = self.dm0
   local ddm0 = self.ddm0
   local m0_new = self.m0_new
   local dm1 = self.dm1
   local ddm1 = self.ddm1
   local m1_new =self.m1_new
   local dm2 = self.dm2
   local ddm2 = self.ddm2
   local m2_new = self.m2_new
   local m2flow = self.m2flow
   local uDrift = self.uDrift
   local vtSq = self.vtSq
   local vtSqTar = self.vtSqTar

   local maxwellian = self.maxwellian
   local calcNumDensity = self.calcNumDensity
   local calcMomDensity = self.calcMomDensity
   local calcKinEnergyDensity = self.calcKinEnergyDensity
   local weakDivide = self.weakDivide
   local weakMultiply = self.weakMultiply
   local weakMultiplyConfPhase = self.weakMultiplyConfPhase

   local vDim = self.vDim

   fOut:copy(fM)

   -- Compute u and vt^2 with the initial M0, M1, M2.
   weakDivide:advance(0., {m0, m1}, {uDrift})
   weakMultiply:advance(0., {uDrift, m1}, {m2flow})
   vtSqTar:combine(1., m2, -1., m2flow)
   weakDivide:advance(0., {m0, vtSqTar}, {vtSqTar})
   vtSqTar:scale(1./vDim)
   local vtSqTar_cur = vtSqTar:get(1)
   vtSqTar:fill(vtSqTar:genIndexer()({2}), vtSqTar_cur)
   local m0_cur = m0:get(1)
   m0:fill(m0:genIndexer()({2}), m0_cur)
   local m1_cur = m1:get(1)
   m1:fill(m1:genIndexer()({2}), m1_cur)
   local m2_cur = m2:get(1)
   m2:fill(m2:genIndexer()({2}), m2_cur)
   --print(m0_cur[1]/math.sqrt(2), m1_cur[1]/math.sqrt(2), m2_cur[1]/math.sqrt(2))
   
   -- Rescale the Maxwellian.
   calcNumDensity:advance(0.0, {fOut}, {m0_new})
   weakDivide:advance(0., {m0_new, m0}, {m0scl})
   weakMultiplyConfPhase:advance(0., {m0scl, fOut}, {fOut})
   -- Compute the moments.
   calcNumDensity:advance(0.0, {fOut}, {m0_new})
   calcMomDensity:advance(0.0, {fOut}, {m1_new})
   calcKinEnergyDensity:advance(0.0, {fOut}, {m2_new})
   
   -- Fix the moments with rescaling and iteration
   local isDone = false
   local step = 0
   while not isDone do
      -- Fix m1.
      ddm1:combine(1., m1, -1., m1_new)
      dm1:combine(1., dm1, 1., ddm1)
      m1_new:combine(1., m1, 1., dm1)
      -- Fix m2.
      ddm2:combine(1., m2, -1., m2_new)
      dm2:combine(1., dm2, 1., ddm2) 
      m2_new:combine(1., m2, 1., dm2)
      local m2_in = m2_new:get(1)
      m2_new:fill(m2_new:genIndexer()({2}), m2_in)
      local m2_in1 = m2_in[1]
 
      -- Compute u and vt^2 with the corrected moments..
      weakDivide:advance(0., {m0_new, m1_new}, {uDrift})
      weakMultiply:advance(0., {uDrift, m1_new}, {m2flow})
      vtSq:combine(1., m2_new, -1., m2flow)
      weakDivide:advance(0., {m0_new, vtSq}, {vtSq})
      vtSq:scale(1./vDim)

      -- Project the Maxwellian onto the basis.
      maxwellian:advance(0.0, {m0_new,uDrift,vtSq,bmag}, {fOut})
      -- Rescale the Maxwellian.
      calcNumDensity:advance(0.0, {fOut}, {m0_new})
      weakDivide:advance(0., {m0_new, m0}, {m0scl})
      weakMultiplyConfPhase:advance(0., {m0scl, fOut}, {fOut})
      -- Compute the moments.
      calcNumDensity:advance(0.0, {fOut}, {m0_new})
      calcMomDensity:advance(0.0, {fOut}, {m1_new})
      calcKinEnergyDensity:advance(0.0, {fOut}, {m2_new})

      local m0_new_cur = m0_new:get(1)
      m0_new:fill(m0_new:genIndexer()({2}), m0_new_cur)
      local m2_new_cur = m2_new:get(1)
      m2_new:fill(m2_new:genIndexer()({2}), m2_new_cur)
      
      -- Compute the l2 norm of the error.
      weakDivide:advance(0., {m0_new, m1_new}, {uDrift})
      weakMultiply:advance(0., {uDrift, m1_new}, {m2flow})
      vtSq:combine(1., m2_new, -1., m2flow)
      weakDivide:advance(0., {m0_new, vtSq}, {vtSq})
      vtSq:scale(1./vDim)
      local err2 = self:l2diff(m2,m2_new) / self:l2norm(m2)
      local err1 = self:l2diff(m1,m1_new) / math.sqrt(self:l2norm(vtSqTar))
      local err = math.max(err1, err2)
      -- Stop at required accuracy.
      if err < self.relEps or step >= self.maxIter then
         isDone = true
         self.flag:appendData(tCurr, {step, err1, err2})
         break
      end
      step = step + 1
      local vtSq_cur = vtSq:get(1)
      --vtSq:fill(vtSq:genIndexer()({2}), vtSq_cur)
      --print(string.format("s=%d, m0_new=%10.8e, m2_in=%10.8e, m2_new=%10.8e, err2=%10.8e", step, m0_new_cur[1]/math.sqrt(2), m2_in1/math.sqrt(2), m2_new_cur[1]/math.sqrt(2), err2))
   end
   t, diagnostic = self.flag:lastData()
end

function IterGkMaxwellianFix:write(tm, frame, speciesName)
   self.flag:write(string.format("%s_%s.bp", speciesName, "convergence"), tm, frame)
end

return IterGkMaxwellianFix
