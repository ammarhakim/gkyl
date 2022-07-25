-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function.
--
-- For LBO collisions this updater also computes the boundary corrections and,
-- if using a piecewise polynomial basis, the star moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local MomDecl      = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local xsys         = require "xsys"

local ffi = require "ffi"
local ffiC = ffi.C
require "Lib.ZeroUtil"

local cudaRunTime
if GKYL_HAVE_CUDA then cudaRunTime = require "Cuda.RunTime" end

ffi.cdef [[ 
// Forward declare for use in function pointers
struct gkyl_mom_type;

/**
 * Function pointer type to compute the needed moment.
 */
typedef void (*momf_t)(const struct gkyl_mom_type *momt,
  const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param);

struct gkyl_mom_type {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of components in moment
  momf_t kernel; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flag;
  struct gkyl_mom_type *on_dev; // pointer to itself or device data
};

/**
 * Create new Vlasov moment type object. Valid 'mom' strings are "M0",
 * "M1i", "M2", "M2ij", "M3i", "M3ijk", "FiveMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

struct gkyl_mom_type* gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

/**
 * Create new Gyrokinetic moment type object. Valid 'mom' strings are "GkM0",
 * "GkM1", "GkM2", "GkM2par", "GkM2perp", "GkM3par", "GkM3perp", "ThreeMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Configuration-space range
 * @param mass Mass of species
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Create new Gyrokinetic moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Set magnitude of the magnetic field, bmag, needed in computing moments.
 * 
 * @param momt Moment type pointer
 * @param bmag Pointer to magnitude of magnetic field
 */
void gkyl_gyrokinetic_set_bmag(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag);

struct gkyl_mom_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_mom_type *momt;

  uint32_t flags;
  struct gkyl_mom_calc *on_dev; // pointer to itself or device data
};

// Object type
typedef struct gkyl_mom_calc gkyl_mom_calc;

/**
 * Create new updater to compute moments of distribution
 * function. Free using gkyl_mom_calc_new_release.
 *
 * @param grid Grid object
 * @param momt Pointer to moment type object
 * @return New updater pointer.
 */
gkyl_mom_calc* gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Create new updater to compute moments of distribution function on
 * NV-GPU. See new() method for documentation.
 */
gkyl_mom_calc* gkyl_mom_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute moment of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Moment calculator updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 * @param mout Output moment array
 */
void gkyl_mom_calc_advance(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);

void gkyl_mom_calc_advance_cu(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);
]]

-- Moments updater object.
local DistFuncMomentCalc = Proto(UpdaterBase)

-- Valid moment names for Vlasov and GK equations.
local goodMomNames = {
   "M0", "M1i", "M2ij", "M2", "M3i", "FiveMoments", "FiveMomentsLBO",
}
local goodGkMomNames = {
   "GkM0", "GkM1", "GkM1proj", "GkM2par", "GkM2perp", "GkM2", "GkM3par", "GkM3perp",
   "GkThreeMoments", "GkThreeMomentsLBO"
}
local goodPartialMomNames = {
   -- Partial velocity moments, integrating over the region
   -- where one of the velocities if positive or negative.
   "M0Pvx",  "M0Pvy",  "M0Pvz", "M0Nvx",  "M0Nvy",  "M0Nvz",
   "M1iPvx", "M1iPvy", "M1iPvz","M1iNvx", "M1iNvy", "M1iNvz",
   "M2Pvx",  "M2Pvy",  "M2Pvz", "M2Nvx",  "M2Nvy",  "M2Nvz",
   "M3iPvx", "M3iPvy", "M3iPvz","M3iNvx", "M3iNvy", "M3iNvz"
}

function DistFuncMomentCalc:isMomentNameGood(nm)
   if lume.find(goodMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isGkMomentNameGood(nm)
   if lume.find(goodGkMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isPartialMomentNameGood(nm)
   if lume.find(goodPartialMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:init(tbl)
   DistFuncMomentCalc.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")
   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   local advArgs = tbl.advanceArgs  -- Sample arguments for advance method.

   self._confBasisID = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Dimension of spaces.
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions.
   self._numBasisP = phaseBasis:numBasis()

   -- Ensure sanity.
   assert(self._polyOrder == phaseBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert((self._confBasisID == phaseBasis:id()) or
          ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and self._confBasisID=="serendipity"),
	  "Type of phase-space and config-space basis must match")

   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   self._fiveMoments, self._fiveMomentsLBO, self.isPartialMom = false, false, false
   if mom == "FiveMoments" or mom == "GkThreeMoments" or
      mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
      self._fiveMoments = true
      if mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
         self._fiveMomentsLBO = true
         -- Rename this variable to call the right kernel below.
         mom = mom=="FiveMomentsLBO" and "FiveMoments" or "GkThreeMoments"
      end
   elseif self:isPartialMomentNameGood(mom) then
      self.isPartialMom = true
   end

   local calcOnDevice = false
   if GKYL_HAVE_CUDA then
      -- This allows us to force an updater to run on the host, even for a GPU simulation.
      self._calcOnHost = tbl.onHost
      if self._calcOnHost then
         self._advanceFunc = self._advance 
      end
   else
      self._calcOnHost = true
   end

   -- Cell index, center and length right of a cell-boundary (also used for current cell for p>1).
   self.idxP = Lin.IntVec(self._pDim)
   self.xcP  = Lin.Vec(self._pDim)
   self.dxP  = Lin.Vec(self._pDim)

   if self.isPartialMom then
      local baseMom
      for _, nm in ipairs(goodMomNames) do 
         baseMom = nm
         if string.find(mom, baseMom) then break end
      end
      -- Extract the direction and whether to integrate over the positive or negative region 
      local velTrans       = {vx=1, vy=2, vz=3}
      local partialMomReg  = string.sub(string.gsub(mom, baseMom, ""),1,1)
      self.partialMomDir   = velTrans[string.sub(string.gsub(mom, baseMom, ""),2)]
      local partialMomDirP = self._cDim+self.partialMomDir
      mom = baseMom
      -- Compute the offsets used to shorten the velocity range. For now assume that the zero
      -- along any velocity dimension is located at a cell boundary and not inside of a cell.
      self.partialMomDirExts   = {0,0}
      local partialMomDirCells = self._onGrid:numCells(partialMomDirP)
      for d = 1,self._pDim do self.idxP[d]=1 end   -- Could be any cell in other directions.
      for idx = 1, partialMomDirCells do
         self.idxP[partialMomDirP] = idx
         self._onGrid:setIndex(self.idxP)
         self._onGrid:cellCenter(self.xcP)
         if (partialMomReg == "P") and (self.xcP[partialMomDirP] > 0.0) then
            self.partialMomDirExts[1] = -(idx-1)
            break
         elseif (partialMomReg == "N") and (self.xcP[partialMomDirP] > 0.0) then
            self.partialMomDirExts[2] = -(partialMomDirCells-(idx-1))
            break
         end
      end
   end

   -- Function to compute specified moment.
   self._isGk = false
   if self:isMomentNameGood(mom) then
      self._kinSpecies = "Vm"
      self._momCalcFun = MomDecl.selectMomCalc(mom, self._confBasisID, self._cDim, self._vDim, self._polyOrder, calcOnDevice)
   elseif self:isGkMomentNameGood(mom) then
      self._kinSpecies = "Gk"
      self._momCalcFun = MomDecl.selectGkMomCalc(mom, self._confBasisID, self._cDim, self._vDim, self._polyOrder)
      self._isGk       = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                           containing the species mass and the background magnetic field
                           to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or FiveMomentsLBO")
   end

   if tbl.gkfacs then
      self.mass    = tbl.gkfacs[1]
      self.bmag    = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
      self.bmagItr = self.bmag:get(1)
   end

   if self._fiveMomentsLBO then
      -- If vDim>1, intFac=2*pi/m or 4*pi/m.
      self._intFac = Lin.Vec(self._vDim)
      for d = 1,self._vDim do self._intFac[d] = 1.0 end
      if self._isGk and (self._vDim > 1) then -- A (vpar,mu) simulation has 3 physical velocity dimensions.
         self._intFac[1] = 2.0*math.pi/self.mass
         self._intFac[2] = 4.0*math.pi/self.mass
      end
      self._isFirst   = true
      self._perpRange = {}    -- Perp ranges in velocity directions.
      if self._polyOrder == 1 then
         self._StarM1iM2Calc = MomDecl.selectStarM1iM2Calc(self._kinSpecies, self._confBasisID, self._cDim, self._vDim)
         -- Cell index, center and length left of a cell-boundary.
         self.idxM = Lin.IntVec(self._pDim)
         self.xcM  = Lin.Vec(self._pDim)
         self.dxM  = Lin.Vec(self._pDim)
      end
   end

   local applyPositivity = xsys.pickBool(tbl.positivity, false)   -- Positivity preserving option.

   self._StarM0Calc     = {} 
   self._uCorrection    = {}
   self._vtSqCorrection = {}
   for vDir = 1, self._vDim do
      if self._isGk and vDir>1 then
         self._StarM0Calc[vDir]  = function(...) end 
         self._uCorrection[vDir] = function(...) end 
      else
         self._StarM0Calc[vDir]  = MomDecl.selectStarM0Calc(vDir, self._kinSpecies, self._confBasisID, self._cDim, self._vDim, self._onGrid:id(), self.applyPositivity)
         self._uCorrection[vDir] = MomDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._confBasisID, self._cDim, self._vDim, self._polyOrder)
      end
      self._vtSqCorrection[vDir] = MomDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._confBasisID, self._cDim, self._vDim, self._polyOrder)
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)

   -- Initialize tools constructed from fields (e.g. ranges).
   self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil

   -- Pre-select and create functions for various cases handled by this updater.
   if self._isGk then 
      self.bmagItrSet = function(cIdx) self.bmag:fill(self.fldTools.confIndexer(cIdx), self.bmagItr) end
   else
      self.bmagItrSet = function(cIdx) end
   end

   if self._fiveMoments then
      -- The first five (Vlasov) or three (GK) moments.
      if self._isGk then   -- Functions that call kernels when computing a single moment.
         self.kernelDriver = function(distfItr, mom0Itr, mom1Itr, mom2Itr)
            self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), 
                             mom0Itr:data(), mom1Itr:data(), mom2Itr:data())
         end
      else
         self.kernelDriver = function(distfItr, mom0Itr, mom1Itr, mom2Itr)
            self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), mom0Itr:data(), mom1Itr:data(), mom2Itr:data())
         end
      end
      if self._fiveMomentsLBO then
         if self._polyOrder == 1 then 
            self.advImpl = function(tCurr, inFld, outFld)
               DistFuncMomentCalc["advanceFiveMomentsLBOp1"](self, tCurr, inFld, outFld)
            end
            if (self._isGk) then
               self.kernelDriverStarM0 = function(vDir, distfItrM, distfItrP, m0StarItr)
                  self._StarM0Calc[vDir](self._intFac[1], self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
               end
               self.kernelDriverStarM1iM2 = function(distfItrP, m1StarItr, m2StarItr)
                  self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), self._intFac[1], self.mass, self.bmagItr:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
               end
            else
               self.kernelDriverStarM0 = function(vDir, distfItrM, distfItrP, m0StarItr)
                  self._StarM0Calc[vDir](self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
               end
               self.kernelDriverStarM1iM2 = function(distfItrP, m1StarItr, m2StarItr)
                  self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
               end
            end
         else
            self.advImpl = function(tCurr, inFld, outFld)
               DistFuncMomentCalc["advanceFiveMomentsLBO"](self, tCurr, inFld, outFld)
            end
         end
         -- Also compute boundary corrections for LBO collisions.
         if self._isGk then   -- Functions that call kernels when computing a single moment.
            self.kernelDriverBoundC = function(vDir, isLo, vBound, distfItrP, cMomBItr, cEnergyBItr)
               self._uCorrection[vDir](isLo, self._intFac[1], vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
               self._vtSqCorrection[vDir](isLo, self._intFac[vDir], vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
            end
         else
            self.kernelDriverBoundC = function(vDir, isLo, vBound, distfItrP, cMomBItr, cEnergyBItr)
               self._uCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
               self._vtSqCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
            end
         end
      else
         self.advImpl = function(tCurr, inFld, outFld)
            DistFuncMomentCalc["advanceFiveMoments"](self, tCurr, inFld, outFld)
         end
      end
   else
      -- Single moment (or partial moment).
      self.advImpl = function(tCurr, inFld, outFld) DistFuncMomentCalc["advanceSingle"](self, tCurr, inFld, outFld) end

      if self._isGk then   -- Functions that call kernels when computing a single moment.
         self.kernelDriver = function(distfItr, momItr)
            self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), momItr:data())
         end
      else
         self.kernelDriver = function(distfItr, momItr)
            self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), momItr:data())
         end
      end
   end

   -- Option to compute moments only once per timestep, based on tCurr input parameter.
   -- NOTE: this should not be used if the updater is used to compute several different quantities in the same timestep.
   self.oncePerTime = xsys.pickBool(tbl.oncePerTime, false)

   if GKYL_USE_GPU then
      if self._isGk then
         self._zero_mom_type = ffiC.gkyl_mom_gyrokinetic_cu_dev_new(confBasis._zero, phaseBasis._zero, self.bmag:localRange(), self.mass, string.sub(mom,3))
         ffiC.gkyl_gyrokinetic_set_bmag(self._zero_mom_type, self.bmag._zeroDevice)
      else
         self._zero_mom_type = ffiC.gkyl_mom_vlasov_cu_dev_new(confBasis._zero, phaseBasis._zero, mom)
      end
      self._zero_mom_calc = ffiC.gkyl_mom_calc_cu_dev_new(self._onGrid._zero, self._zero_mom_type)
   else
      if self._isGk then
         self._zero_mom_type = ffiC.gkyl_mom_gyrokinetic_new(confBasis._zero, phaseBasis._zero, self.bmag:localRange(), self.mass, string.sub(mom,3))
         ffiC.gkyl_gyrokinetic_set_bmag(self._zero_mom_type, self.bmag._zero)
      else
         self._zero_mom_type = ffiC.gkyl_mom_vlasov_new(confBasis._zero, phaseBasis._zero, mom)
      end
      self._zero_mom_calc = ffiC.gkyl_mom_calc_new(self._onGrid._zero, self._zero_mom_type)
   end
end

function DistFuncMomentCalc:initFldTools(inFld, outFld)
   -- Pre-initialize tools (ranges, pointers, etc) depending on fields and used in the advance method.
   local tools = {}

   local distf, mom = inFld[1], outFld[1]

   local cDim, vDim = self._cDim, self._vDim

   local phaseRange = distf:localRange()
   if self.onGhosts then -- Extend range to config-space ghosts.
      for dir = 1, cDim do 
         phaseRange = phaseRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
      end
   end

   -- Construct ranges for nested loops.
   tools.confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = self._onGrid:numSharedProcs() }
   tools.velRange = phaseRange:selectLast(vDim)
   tools.tId      = self._onGrid:subGridSharedId()    -- Local thread ID.

   if self.isPartialMom then
      -- For partial moments, we only wish to integrate over part of one of the velocity dimensions.
      tools.velRange = tools.velRange:extendDir(self.partialMomDir,self.partialMomDirExts[1],self.partialMomDirExts[2])
   end

   tools.phaseIndexer, tools.confIndexer = distf:genIndexer(), mom:genIndexer()

   if self._fiveMomentsLBO then
      tools.perpRange, tools.dirLoIdx, tools.dirUpIdx = {}, {}, {}
      for vDir = 1, vDim do
         -- Lower/upper bounds in direction 'vDir': cell indices.
         if self._polyOrder==1 then
            -- Including outer edges.
            tools.dirLoIdx[vDir], tools.dirUpIdx[vDir] = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)+1
         else
            tools.dirLoIdx[vDir], tools.dirUpIdx[vDir] = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)
         end

         tools.perpRange[vDir] = phaseRange
         -- Shorten configuration range.
         for cd = 1, cDim do tools.perpRange[vDir] = tools.perpRange[vDir]:shorten(cd) end
         tools.perpRange[vDir] = tools.perpRange[vDir]:shorten(cDim+vDir) -- Velocity range orthogonal to 'vDir'.
      end
   end

   return tools
end

function DistFuncMomentCalc:advanceSingle(tCurr, inFld, outFld)
   -- Compute a single moment of the distribution function.
   local distf, mom = inFld[1], outFld[1]

   local distfItr, momItr = distf:get(1), mom:get(1)
   
   mom._zero:clear(0.0) -- Zero out moments.

   local grid = self._onGrid

   -- Outer loop is threaded and over configuration space.
   for cIdx in self.fldTools.confRangeDecomp:rowMajorIter(self.fldTools.tId) do
   
      cIdx:copyInto(self.idxP)

      self.bmagItrSet(cIdx)
      mom:fill(self.fldTools.confIndexer(cIdx), momItr)
   
      -- Inner loop is over velocity space: no threading to avoid race conditions.
      for vIdx in self.fldTools.velRange:rowMajorIter() do
         for d = 1, self._vDim do self.idxP[self._cDim+d] = vIdx[d] end
      
         grid:setIndex(self.idxP)
         grid:cellCenter(self.xcP)
         grid:getDx(self.dxP)
      
         distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItr)

         self.kernelDriver(distfItr, momItr)
      end
   end
end

function DistFuncMomentCalc:advanceFiveMoments(tCurr, inFld, outFld)
   -- Compute the first five (Vlasov) or three (GK) moments of the distribution function.
   local distf            = inFld[1]
   local mom0, mom1, mom2 = outFld[1], outFld[2], outFld[3] 

   local distfItr                  = distf:get(1)
   local mom0Itr, mom1Itr, mom2Itr = mom0:get(1), mom1:get(1), mom2:get(1) 
   
   mom0._zero:clear(0.0);  mom1._zero:clear(0.0);  mom2._zero:clear(0.0) -- Zero out moments.
   
   local grid = self._onGrid

   -- Outer loop is threaded and over configuration space.
   for cIdx in self.fldTools.confRangeDecomp:rowMajorIter(self.fldTools.tId) do
   
      cIdx:copyInto(self.idxP)

      self.bmagItrSet(cIdx)
      mom0:fill(self.fldTools.confIndexer(cIdx), mom0Itr)
      mom1:fill(self.fldTools.confIndexer(cIdx), mom1Itr)
      mom2:fill(self.fldTools.confIndexer(cIdx), mom2Itr)
   
      -- Inner loop is over velocity space: no threading to avoid race conditions.
      for vIdx in self.fldTools.velRange:rowMajorIter() do
         for d = 1, self._vDim do self.idxP[self._cDim+d] = vIdx[d] end
      
         grid:setIndex(self.idxP)
         grid:cellCenter(self.xcP)
         grid:getDx(self.dxP)
      
         distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItr)

         self.kernelDriver(distfItr, mom0Itr, mom1Itr, mom2Itr)
      end
   end
end

function DistFuncMomentCalc:advanceFiveMomentsLBO(tCurr, inFld, outFld)
   -- Compute the first five (Vlasov) or three (GK) moments of the distribution function,
   -- and also calculate the boundary corrections for LBO collisions.
   local distf            = inFld[1]
   local mom0, mom1, mom2 = outFld[1], outFld[2], outFld[3]
   local cMomB, cEnergyB  = outFld[4], outFld[5] 

   local distfItr                  = distf:get(1)
   local mom0Itr, mom1Itr, mom2Itr = mom0:get(1), mom1:get(1), mom2:get(1)
   local cMomBItr, cEnergyBItr     = cMomB:get(1), cEnergyB:get(1) 
   -- For corrections and star moments: distribution functions left/right of a cell-boundary.
   local distfItrP, distfItrM      = distf:get(1), distf:get(1)
   
   mom0._zero:clear(0.0);  mom1._zero:clear(0.0);  mom2._zero:clear(0.0) -- Zero out moments.
   cMomB._zero:clear(0.0);  cEnergyB._zero:clear(0.0)
   
   local grid       = self._onGrid
   local cDim, vDim = self._cDim, self._vDim

   -- Outer loop is threaded and over configuration space.
   for cIdx in self.fldTools.confRangeDecomp:rowMajorIter(self.fldTools.tId) do
   
      cIdx:copyInto(self.idxP)

      self.bmagItrSet(cIdx)
      mom0:fill(self.fldTools.confIndexer(cIdx), mom0Itr)
      mom1:fill(self.fldTools.confIndexer(cIdx), mom1Itr)
      mom2:fill(self.fldTools.confIndexer(cIdx), mom2Itr)
   
      -- Inner loop is over velocity space: no threading to avoid race conditions.
      for vIdx in self.fldTools.velRange:rowMajorIter() do
         for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end
      
         grid:setIndex(self.idxP)
         grid:cellCenter(self.xcP)
         grid:getDx(self.dxP)
      
         distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItr)
   
         self.kernelDriver(distfItr, mom0Itr, mom1Itr, mom2Itr)
      end
   
      -- Now loop over velocity space boundary surfaces to compute boundary corrections.
      cMomB:fill(self.fldTools.confIndexer(cIdx), cMomBItr)
      cEnergyB:fill(self.fldTools.confIndexer(cIdx), cEnergyBItr)
   
      -- isLo=true current cell is the lower boundary cell.
      -- isLo=false current cell is the upper boundary cell.
      local isLo = true
   
      for vDir = 1, vDim do
         for vPerpIdx in self.fldTools.perpRange[vDir]:rowMajorIter() do
            vPerpIdx:copyInto(self.idxP)
            for d = 1, cDim do self.idxP[d] = cIdx[d] end
   
            for _, i in ipairs({self.fldTools.dirLoIdx[vDir], self.fldTools.dirUpIdx[vDir]}) do   -- This loop is over edges.
               self.idxP[cDim+vDir] = i
   
               grid:setIndex(self.idxP)
               grid:getDx(self.dxP)
               grid:cellCenter(self.xcP)
   
               distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItrP)
   
               local vBound = isLo and grid:cellLowerInDir(cDim + vDir) or grid:cellUpperInDir(cDim + vDir)
   
               self.kernelDriverBoundC(vDir, isLo, vBound, distfItrP, cMomBItr, cEnergyBItr)
   
               isLo = not isLo
            end
         end    -- vPerpIdx loop.
      end    -- vDir loop.

   end    -- Loop over configuration space.
end

function DistFuncMomentCalc:advanceFiveMomentsLBOp1(tCurr, inFld, outFld)
   -- Like advanceFiveMomentsLBO but for polyOrder=1, i.e. it also computes the star moments.
   local distf                  = inFld[1]
   local mom0, mom1, mom2       = outFld[1], outFld[2], outFld[3]
   local cMomB, cEnergyB        = outFld[4], outFld[5] 
   local m0Star, m1Star, m2Star = outFld[6], outFld[7], outFld[8]

   local distfItrP, distfItrM            = distf:get(1), distf:get(1)
   local mom0Itr, mom1Itr, mom2Itr       = mom0:get(1), mom1:get(1), mom2:get(1)
   local cMomBItr, cEnergyBItr           = cMomB:get(1), cEnergyB:get(1) 
   local m0StarItr, m1StarItr, m2StarItr = m0Star:get(1), m1Star:get(1), m2Star:get(1)  
   
   mom0._zero:clear(0.0);  mom1._zero:clear(0.0);  mom2._zero:clear(0.0) -- Zero out moments.
   cMomB._zero:clear(0.0);  cEnergyB._zero:clear(0.0)
   m0Star._zero:clear(0.0);  m1Star._zero:clear(0.0);  m2Star._zero:clear(0.0)
   
   local grid       = self._onGrid
   local cDim, vDim = self._cDim, self._vDim

   -- Outer loop is threaded and over configuration space.
   for cIdx in self.fldTools.confRangeDecomp:rowMajorIter(self.fldTools.tId) do
   
      cIdx:copyInto(self.idxP)
   
      self.bmagItrSet(cIdx)
      mom0:fill(self.fldTools.confIndexer(cIdx), mom0Itr)
      mom1:fill(self.fldTools.confIndexer(cIdx), mom1Itr)
      mom2:fill(self.fldTools.confIndexer(cIdx), mom2Itr)
      cMomB:fill(self.fldTools.confIndexer(cIdx), cMomBItr)
      cEnergyB:fill(self.fldTools.confIndexer(cIdx), cEnergyBItr)
      -- To have energy conservation with piece-wise linear, we must use star moments
      -- in the second equation of the weak system solved in SelfPrimMoments.
      m0Star:fill(self.fldTools.confIndexer(cIdx), m0StarItr)
      m1Star:fill(self.fldTools.confIndexer(cIdx), m1StarItr)
      m2Star:fill(self.fldTools.confIndexer(cIdx), m2StarItr)
   
      -- Only when the contributions to m0Star from the first direction are collected, do we
      -- collect contributions to m1Star and m2Star. Also, since Gk velocities are organized as
      -- (vpar,mu) the velocity correction is only computed for the first velocity direction.
      local firstDir = true
   
      local isLo = true   -- isLo=true/false current cell is the lower/upper boundary cell.
   
      for vDir = 1, vDim do
         -- Outer loop is over directions orthogonal to 'vDir' and
         -- inner loop is over 1D slice in 'vDir'.
         for vPerpIdx in self.fldTools.perpRange[vDir]:rowMajorIter() do
            vPerpIdx:copyInto(self.idxM); vPerpIdx:copyInto(self.idxP)
            for d = 1, cDim do 
               self.idxM[d] = cIdx[d];  self.idxP[d] = cIdx[d]
            end
   
            for i = self.fldTools.dirLoIdx[vDir], self.fldTools.dirUpIdx[vDir] do     -- This loop is over edges.
               self.idxM[cDim+vDir], self.idxP[cDim+vDir] = i-1, i -- Cell left/right of edge 'i'.
   
               grid:setIndex(self.idxM);  grid:getDx(self.dxM);  grid:cellCenter(self.xcM);
               grid:setIndex(self.idxP);  grid:getDx(self.dxP);  grid:cellCenter(self.xcP)

               distf:fill(self.fldTools.phaseIndexer(self.idxM), distfItrM)
               distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItrP)
   
               if i>self.fldTools.dirLoIdx[vDir] and i<self.fldTools.dirUpIdx[vDir] then
                  self.kernelDriverStarM0(vDir, distfItrM, distfItrP, m0StarItr)
               end

               if firstDir and i<self.fldTools.dirUpIdx[vDir] then
                  self.kernelDriver(distfItrP, mom0Itr, mom1Itr, mom2Itr)
                  self.kernelDriverStarM1iM2(distfItrP, m1StarItr, m2StarItr)
               end
   
               if i==self.fldTools.dirLoIdx[vDir] or i==self.fldTools.dirUpIdx[vDir]-1 then
                  -- CAREFUL: for vBound below we assume idxP was set after idxM above.
                  local vBound = isLo and grid:cellLowerInDir(cDim + vDir) or grid:cellUpperInDir(cDim + vDir)

                  self.kernelDriverBoundC(vDir, isLo, vBound, distfItrP, cMomBItr, cEnergyBItr)
   
                  isLo = not isLo
               end
   
            end    -- Loop over edges.
         end    -- Loop over directions perpendicular to vDir.
         firstDir = false
   
      end    -- vDir loop.
   end    -- Loop over configuration space.
end

-- Advance method.
function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   if self._zero_mom_calc then
      if self.oncePerTime and self.tCurr == tCurr then return end -- Do nothing, already computed on this step.
      local distf = inFld[1]
      local mout = outFld[1]
      local phaseRange = distf:localRange()
      local confRange
      -- need to use localExtRange for confRange when onGhosts=true
      -- note that phaseRange does not need to use localExtRange because
      -- phaseRange is only used to get the velocity-space sub-range,
      -- which should not include ghosts
      if self.onGhosts then
         confRange = mout:localExtRange()
      else
         confRange = mout:localRange()
      end
      ffiC.gkyl_mom_calc_advance(self._zero_mom_calc, phaseRange, confRange, distf._zero, mout._zero)
      if self.oncePerTime then self.tCurr = tCurr end
   else
      if self.oncePerTime and self.tCurr == tCurr then return end -- Do nothing, already computed on this step.

      self.fldTools = self.fldTools or self:initFldTools(inFld,outFld)

      self.advImpl(tCurr, inFld, outFld)

      if self.oncePerTime then self.tCurr = tCurr end
   end
end

function DistFuncMomentCalc:_advanceOnDevice(tCurr, inFld, outFld)
   if self.oncePerTime and self.tCurr == tCurr then return end -- Do nothing, already computed on this step.
   local distf = inFld[1]
   local mout = outFld[1]
   local phaseRange = distf:localRange()
   local confRange
   if self.onGhosts then
      -- on GPU, we do need to set up phaseRange with ghosts in conf 
      -- space but not in v-space as a subrange of phaseExtRange
      local phaseExtRange = distf:localExtRange()
      local lower = phaseExtRange:lowerAsVec()
      local upper = phaseExtRange:upperAsVec()
      for d=self._cDim+1, self._pDim do
         lower[d] = lower[d] + distf:lowerGhost()
         upper[d] = upper[d] - distf:upperGhost()
      end
      phaseRange = phaseExtRange:subRange(lower, upper)
      confRange = mout:localExtRange()
   else
      confRange = mout:localRange()
   end
   mout:clear(0.0)
   ffiC.gkyl_mom_calc_advance_cu(self._zero_mom_calc, phaseRange, confRange, distf._zeroDevice, mout._zeroDevice)
   if self.oncePerTime then self.tCurr = tCurr end
end


return DistFuncMomentCalc
