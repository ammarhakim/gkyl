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
local IterMaxwellianFix = Proto(UpdaterBase)

-- constructor
function IterMaxwellianFix:init(tbl)
   IterMaxwellianFix.super.init(self, tbl) -- setup base object

   self.phaseGrid =
      assert(tbl.onGrid, "Updater.IterMaxwellianFix: Must provide phase space grid object 'onGrid'")
   self.confGrid =
      assert(tbl.confGrid, "Updater.IterMaxwellianFix: Must provide configuration space grid object 'confGrid'")
   self.confBasis =
      assert(tbl.confBasis, "Updater.IterMaxwellianFix: Must provide configuration space basis object 'confBasis'")
   self.phaseBasis =
      assert(tbl.phaseBasis, "Updater.IterMaxwellianFix: Must provide phase space basis object 'phaseBasis'")

   self.maxIter = tbl.maxIter
   self.relEps = tbl.relEps
   
   local phaseGrid, confGrid, confBasis, phaseBasis =
      self.phaseGrid, self.confGrid, self.confBasis, self.phaseBasis

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

   -- Allocate memory
   self.m0scl = getField(confGrid, confBasis, 1) 
   self.dm0 = getField(confGrid, confBasis, 1) 
   self.ddm0 = getField(confGrid, confBasis, 1) 
   self.m0_new = getField(confGrid, confBasis, 1)
   self.dm1i = getField(confGrid, confBasis, self.vDim) 
   self.ddm1i = getField(confGrid, confBasis, self.vDim) 
   self.m1i_new = getField(confGrid, confBasis, self.vDim)
   self.dm2 = getField(confGrid, confBasis, 1)
   self.ddm2 = getField(confGrid, confBasis, 1)
   self.m2_new = getField(confGrid, confBasis, 1)
   self.m2flow = getField(confGrid, confBasis, 1)
   self.uDrift = getField(confGrid, confBasis, self.vDim)
   self.vtSq = getField(confGrid, confBasis, 1)

   -- create updaters
   self.maxwellian = MaxwellianOnBasis {
      onGrid = phaseGrid,
      confGrid = confGrid,
      phaseBasis = phaseBasis,
      confBasis = confBasis,
   }

   -- Moment updaters.
   self.calcNumDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "M0",
   }
   self.calcMomDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "M1i",
   }
   self.calcKinEnergyDensity = DistFuncMomentCalc {
      onGrid = phaseGrid,
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      moment = "M2",
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
end

-- compute the L2 norm
function IterMaxwellianFix:l2norm(fld)
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

function IterMaxwellianFix:l2diff(f1, f2)
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

-- advance method
function IterMaxwellianFix:_advance(tCurr, inFld, outFld)

   -- fM is the uncorrected Maxwellian projection
   local fM = assert(inFld[1], "IterMaxwellianFix.advance: Must-specify an input distribution function")
   local m0 = assert(inFld[2], "IterMaxwellianFix.advance: Must-specify an input number density")
   local m1i = assert(inFld[3], "IterMaxwellianFix.advance: Must-specify an input moment density")
   local m2 = assert(inFld[4], "IterMaxwellianFix.advance: Must-specify an input energy density")
   
   local fOut = assert(outFld[1], "IterMaxwellianFix.advance: Must-specify an output distribution function")

   local m0scl = self.m0scl
   local dm0 = self.dm0
   local ddm0 = self.ddm0
   local m0_new = self.m0_new
   local dm1i = self.dm1i
   local ddm1i = self.ddm1i
   local m1i_new =self.m1i_new
   local dm2 = self.dm2
   local ddm2 = self.ddm2
   local m2_new = self.m2_new
   local m2flow = self.m2flow
   local uDrift = self.uDrift
   local vtSq = self.vtSq

   local maxwellian = self.maxwellian
   local calcNumDensity = self.calcNumDensity
   local calcMomDensity = self.calcMomDensity
   local calcKinEnergyDensity = self.calcKinEnergyDensity
   local weakDivide = self.weakDivide
   local weakMultiply = self.weakMultiply
   local weakMultiplyConfPhase = self.weakMultiplyConfPhase

   local vDim = self.vDim

   fOut:copy(fM)

   -- Rescale the Maxwellian.
   calcNumDensity:advance(0.0, {fOut}, {m0_new})
   weakDivide:advance(0., {m0_new, m0}, {m0scl})
   weakMultiplyConfPhase:advance(0., {m0scl, fOut}, {fOut})
   -- Compute the moments.
   calcNumDensity:advance(0.0, {fOut}, {m0_new})
   calcMomDensity:advance(0.0, {fOut}, {m1i_new})
   calcKinEnergyDensity:advance(0.0, {fOut}, {m2_new})
   
   -- Fix the moments with rescaling and iteration
   local isDone = false
   local step = 0
   while not isDone do
      -- Fix m1i.
      ddm1i:combine(1., m1i, -1., m1i_new)
      dm1i:combine(1., dm1i, 1., ddm1i)
      m1i_new:combine(1., m1i, 1., dm1i)
      -- Fix m2.
      ddm2:combine(1., m2, -1., m2_new)
      dm2:combine(1., dm2, 1., ddm2) 
      m2_new:combine(1., m2, 1., dm2)
      
      -- Compute u and vt^2 with the moments of the rescaled Maxwellian.
      weakDivide:advance(0., {m0_new, m1i_new}, {uDrift})
      weakMultiply:advance(0., {uDrift, m1i_new}, {m2flow})
      vtSq:combine(1., m2_new, -1., m2flow)
      weakDivide:advance(0., {m0_new, vtSq}, {vtSq})
      vtSq:scale(1./vDim)
      
      -- Project the Maxwellian onto the basis.
      maxwellian:advance(0.0, {m0_new,uDrift,vtSq}, {fOut})
      -- Rescale the Maxwellian.
      calcNumDensity:advance(0.0, {fOut}, {m0_new})
      weakDivide:advance(0., {m0_new, m0}, {m0scl})
      weakMultiplyConfPhase:advance(0., {m0scl, fOut}, {fOut})
      -- Compute the moments.
      calcNumDensity:advance(0.0, {fOut}, {m0_new})
      calcMomDensity:advance(0.0, {fOut}, {m1i_new})
      calcKinEnergyDensity:advance(0.0, {fOut}, {m2_new})
      
      -- Compute the l2 norm of the truncation error.
      local err = self:l2diff(m2,m2_new) / self:l2norm(m2)
      -- Stop at required accuracy.
      if err < self.relEps or step >= self.maxIter then
         isDone = true
         break
      end
      step = step + 1
      print("Step number:", step)
   end
   
end

return IterMaxwellianFix
