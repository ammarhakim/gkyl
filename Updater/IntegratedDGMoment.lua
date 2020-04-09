-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated moments of a DG field.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc        = require "Lib.Alloc"
local Lin          = require "Lib.Linalg"
local IntMomDecl   = require "Updater.momentCalcData.IntegratedDGMomentModDecl"
local Mpi          = require "Comm.Mpi"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"

-- Integrated moments updater object.
local IntegratedDGMoment = Proto(UpdaterBase)

local goodMomNames = {
   -- In the descriptions below int() means integrate over all space. 
   -- As needed, other developers/users may add the functionality to compute multiple
   -- of these at once, so that the grid doesn't need to be traversed multiple times.
   "one",                                            -- int(f)  
   "x1", "x2", "x3", "x4", "x5", "x6",               -- int(x_i*f)
   "x1Sq", "x2Sq", "x3Sq", "x4Sq", "x5Sq", "x6Sq",   -- int(x_i^2*f)
   "xi",                                             -- int(x_i*f) for i=1,2,...,dim.
   "xSq",                                            -- int((sum_i x_i^2)*f)
   "v1", "v2", "v3",                                 -- int(v_i*f), only for phase space fields. 
   "v1Sq", "v2Sq", "v3Sq",                           -- int(v_i^2*f), only for phase space fields.  
   "vi",                                             -- Compute "v1", "v2", "v3" all at once.
   "vSq",                                            -- int((sum_i v_i^2)*f), only for phase space fields.  
   "intM"                                            -- Compute "one", "v1", "v2", "v3", "vSq", all at once.
}

local goodVspaceMomNames = {
   -- In the descriptions below int() means integrate over all space. 
   "v1", "v2", "v3",                                 -- int(v_i*f), only for phase space fields. 
   "v1Sq", "v2Sq", "v3Sq",                           -- int(v_i^2*f), only for phase space fields.  
   "vSq",                                            -- int((sum_i v_i^2)*f), only for phase space fields.  
   "vi",                                             -- Compute "v1", "v2", "v3" all at once.
   "intM"                                            -- Compute "one", "v1", "v2", "v3", "vSq", all at once.
}

local function isMomentNameGood(nm)
   if lume.find(goodMomNames, nm) then
      return true
   end
   return false
end

local function isVspaceMoment(nm)
   if lume.find(goodVspaceMomNames, nm) then
      return true
   end
   return false
end

function IntegratedDGMoment:init(tbl)
   IntegratedDGMoment.super.init(self, tbl) -- Setup base object.

   self.onGrid   = assert(
      tbl.onGrid, "Updater.IntegratedDGMoment: Must provide grid object using 'onGrid'")
   local onBasis = assert(
      tbl.basis,
      "Updater.IntegratedDGMoment: Must provide basis object using 'onBasis'")
   local mom     = assert(
      tbl.moment, "Updater.IntegratedDGMoment: Must provide moment to compute using 'moment'.")

   local basisID, polyOrder = onBasis:id(), onBasis:polyOrder()

   if (not isMomentNameGood(mom)) then
      assert(false, "Updater.IntegratedDGMoment: invalid moment given in 'moment'. See valid list in updater.")
   end

   self.dim = onBasis:ndim()

   self.dxv = Lin.Vec(self.dim)     -- Cell shape.
   self.w   = Lin.Vec(self.dim)     -- Cell center.

   local momDirIdxS = string.find(mom,'%d')
   local momDir     = nil
   if (momDirIdxS) then momDir = tonumber(string.sub(mom,momDirIdxS,momDirIdxS)) end


   self.intMomNcomp = 1      -- Number of integrated moments computed.

   if isVspaceMoment(mom) then
      local confBasis = assert(tbl.confBasis,
                               "Updater.IntegratedDGMoment: Must provide configuration-space basis object using 'confBasis'")
      -- Ensure sanity.
      assert(polyOrder == confBasis:polyOrder(),
             "Updater.IntegratedDGMoment: Polynomial orders of phase-space and config-space basis must match")
      assert(basisID == confBasis:id(),
             "Updater.IntegratedDGMoment: Type of phase-space and config-space basis must match")

      -- Dimension of configuration and velocity spaces.
      self.cDim = confBasis:ndim()
      self.vDim = self.dim - self.cDim

      if (momDirIdxS) then
         if (momDir < 1) or (momDir > self.vDim) then
            assert(false, "Updater.IntegratedDGMoment: velocity coordinate given in 'moment' is invalid for this dimensionality. See valid list in updater.")
         end
      end

     if (mom == 'vi') then
        -- This computes the integrals of v_1, v_2 and v_3 at the same time.
        self.intMomNcomp = self.vDim
     elseif (mom == 'intM') then
        -- This computes the integrals of M0, M1i and M2 at the same time.
        self.intMomNcomp = 2+self.vDim
     end

      -- Kernel to compute integrated velocity moment.
      self.intMomCalc = IntMomDecl.selectIntVspaceMomCalc(mom, basisID, self.cDim, self.vDim, polyOrder)

   else
      if (momDirIdxS) then
         if (momDir < 1) or (momDir > self.dim) then
            assert(false, "Updater.IntegratedDGMoment: coordinate given in 'moment' is invalid for this dimensionality. See valid list in updater.")
         end
      end

     if (mom == 'xi') then
        -- This computes the integrals of x_1, x_2, etc, at the same time.
        self.intMomNcomp = self.dim
     end

      -- Kernel to compute integrated moment.
      self.intMomCalc = IntMomDecl.selectIntMomCalc(mom, basisID, self.dim, polyOrder)

   end

   self.localMom  = Lin.Vec(self.intMomNcomp)
   self.globalMom = Lin.Vec(self.intMomNcomp)

end

-- Advance method.
function IntegratedDGMoment:_advance(tCurr, inFld, outFld)
   local grid           = self.onGrid
   local scalarFld, mom = inFld[1], outFld[1]

   assert(mom:numComponents() == self.intMomNcomp,
          "Updater.IntegratedDGMoment: numComponents in DynVector is wrong for the moment requested.")

   local fldIndexer = scalarFld:genIndexer()
   local fldItr     = scalarFld:get(1)

   -- Clear local values.
   for i = 1,self.intMomNcomp do self.localMom[i] = 0.0 end

   -- Construct range for shared memory.
   local fldRange       = scalarFld:localRange()
   local fldRangeDecomp = LinearDecomp.LinearDecompRange {
      range = fldRange:selectFirst(self.dim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.

   -- Loop, computing integrated moments in each cell.
   for idx in fldRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(idx)
      grid:cellCenter(self.w)
      for d = 1, self.dim do self.dxv[d] = grid:dx(d) end
      scalarFld:fill(fldIndexer(idx), fldItr)
      self.intMomCalc(self.w:data(), self.dxv:data(), fldItr:data(), self.localMom:data())
   end

   -- All-reduce across processors and push result into dyn-vector.
   Mpi.Allreduce(
      self.localMom:data(), self.globalMom:data(), self.intMomNcomp, Mpi.DOUBLE, Mpi.SUM, self:getComm())

   mom:appendData(tCurr, self.globalMom)
end

return IntegratedDGMoment
