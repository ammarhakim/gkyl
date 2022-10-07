-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the Spitzer collision frequency.
-- Eventually this updater may:
-- a) Return a cell-wise constant nu given a normalized initial normNu.
-- b) Return a nu expanded in the basis function given an initial normNu.
-- c) Build nu in SI units (either cell-wise constant or varying).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase   = require "Updater.Base"
local Proto         = require "Lib.Proto"
local SpitzerNuDecl = require "Updater.spitzerNuCalcData.SpitzerNuModDecl"
local xsys          = require "xsys"
local Lin           = require "Lib.Linalg"
local ffi           = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_spitzer_coll_freq gkyl_spitzer_coll_freq;

/**
 * Create new updater to either compute the Spitzer collision frequency from
 * scratch based on local parameters, or scale a normalized collision frequency
 * by the local n_r/(v_ts^2+v_tr^2)^(3/2).
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_spitzer_coll_freq* gkyl_spitzer_coll_freq_new(
  const struct gkyl_basis *basis, int num_quad, bool use_gpu);

/**
 * Scale the normalized collision frequency, normNu, by
 * n_r/(v_ts^2+v_tr^2)^(3/2) and project it on to the basis.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range
 * @param vtSqSelf Thermal speed squared of this species.
 * @param m0Other Thermal speed squared of the other species.
 * @param vtSqOther Thermal speed squared of the other species.
 * @param normNu Normalized collision frequency to scale.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance_normnu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  double normNu, struct gkyl_array *nuOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_spitzer_coll_freq_release(gkyl_spitzer_coll_freq* up);
]]

-- Updater object.
local SpitzerCollisionality = Proto(UpdaterBase)

function SpitzerCollisionality:init(tbl)
   SpitzerCollisionality.super.init(self, tbl) -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SpitzerCollisionality: Must provide grid object using 'onGrid'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SpitzerCollisionality: Must provide the configuration basis object using 'confBasis'.")

   self._cellConstNu = assert(
      tbl.useCellAverageNu, "Updater.SpitzerCollisionality: Must specify whether to use cell averaged collisionality in 'useCellAverageNu'.")

   self._isNormNu = tbl.willInputNormNu
   assert(self._isNormNu~=nil, "Updater.SpitzerCollisionality: Must specify whether normalized collisionality (normNu) will be inputed or not.")

   self._elemCharge = assert(
      tbl.elemCharge, "Updater.SpitzerCollisionality: Must specify elementary charge ('elemCharge') to build Spitzer collisionality.")
   self._epsilon0   = assert(
      tbl.epsilon0, "Updater.SpitzerCollisionality: Must specify vacuum permittivity ('epsilon0') to build Spitzer collisionality.")
   self._hBar       = assert(
      tbl.hBar, "Updater.SpitzerCollisionality: Must specify Planck's constant h divided by 2pi ('hBar') to build Spitzer collisionality.")
   local nuFracIn = tbl.nuFrac
   if nuFracIn then
      self._nuFrac = nuFracIn
   else
      self._nuFrac = 1.0
   end

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   if self._isNormNu then
      local numQuad = confBasis:polyOrder()+1
      self._zero = ffi.gc(ffiC.gkyl_spitzer_coll_freq_new(confBasis._zero, numQuad, self._useGPU),
                          ffiC.gkyl_spitzer_coll_freq_release)
      self.advanceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advance_normNu"](self, tCurr, inFld, outFld)
      end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advanceOnDevice_normNu"](self, tCurr, inFld, outFld)
      end
   else

      -- Dimension of configuration space.
      local cDim = confBasis:ndim()
      -- Basis name and polynomial order.
      local basisID   = confBasis:id()
      local polyOrder = confBasis:polyOrder()
      -- Number of basis functions.
      local numBasisC = confBasis:numBasis()

      if self._isNormNu then
         if self._cellConstNu then
            self._SpitzerNuCalc = SpitzerNuDecl.selectCellAvSpitzerNuScale(basisID, cDim, polyOrder)
         else
            self._SpitzerNuCalc = SpitzerNuDecl.selectSpitzerNuScale(basisID, cDim, polyOrder)
         end
      else
         if self._cellConstNu then
            self._SpitzerNuCalc = SpitzerNuDecl.selectCellAvSpitzerNuBuild(basisID, cDim, polyOrder)
         else
            self._SpitzerNuCalc = SpitzerNuDecl.selectSpitzerNuBuild(basisID, cDim, polyOrder)
         end
      end
      self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

      self._BmagZero = Lin.Vec(numBasisC)
      for iC = 1,numBasisC do self._BmagZero[iC] = 0.0 end
   
      self.advanceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advance_build"](self, tCurr, inFld, outFld)
      end
   end

end

-- Advance methods.
function SpitzerCollisionality:_advance_normNu(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local normNu = inFld[9]
   local Bmag   = inFld[10]

   local nuOut = outFld[1]

   local range = self.onGhost and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance_normnu(self._zero, range,
      vtSqSelf._zero, m0Other._zero, vtSqOther._zero, normNu*self._nuFrac, nuOut._zero)
end

function SpitzerCollisionality:_advanceOnDevice_normNu(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local normNu = inFld[9]
   local Bmag   = inFld[10]

   local nuOut = outFld[1]

   local range = self.onGhost and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance_normnu(self._zero, range,
      vtSqSelf._zeroDevice, m0Other._zeroDevice, vtSqOther._zeroDevice, normNu*self._nuFrac, nuOut._zeroDevice)
end

function SpitzerCollisionality:_advance_build(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local normNu = inFld[9]
   local Bmag   = inFld[10]

   local nuOut = outFld[1]

   local confIndexer = m0Self:genIndexer()

   local m0SelfItr    = m0Self:get(1)
   local vtSqSelfItr  = vtSqSelf:get(1)
   local m0OtherItr   = m0Other:get(1)
   local vtSqOtherItr = vtSqOther:get(1)

   local nuOutItr = nuOut:get(1)

   local confRange = m0Self:localRange()
   if self.onGhosts then confRange = m0Self:localExtRange() end

   local grid = self._onGrid

   -- Fork logic here to avoid if-statement in space loop.
   if (Bmag ~= nil) then
      local BmagItr = Bmag:get(1)
      for cIdx in confRange:rowMajorIter() do
         grid:setIndex(cIdx)
 
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)

         Bmag:fill(confIndexer(cIdx), BmagItr)
 
         nuOut:fill(confIndexer(cIdx), nuOutItr)
 
         self._SpitzerNuCalc(self._elemCharge, self._epsilon0, self._hBar, self._nuFrac, chargeSelf, massSelf, m0SelfItr:data(), vtSqSelfItr:data(), chargeOther, massOther, m0OtherItr:data(), vtSqOtherItr:data(), normNu, BmagItr:data(), nuOutItr:data())
      end
   else
      local BmagItr = self._BmagZero
      for cIdx in confRange:rowMajorIter() do
         grid:setIndex(cIdx)
 
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)
 
         nuOut:fill(confIndexer(cIdx), nuOutItr)
 
         self._SpitzerNuCalc(self._elemCharge, self._epsilon0, self._hBar, self._nuFrac, chargeSelf, massSelf, m0SelfItr:data(), vtSqSelfItr:data(), chargeOther, massOther, m0OtherItr:data(), vtSqOtherItr:data(), normNu, BmagItr:data(), nuOutItr:data())
      end
   end
end

function SpitzerCollisionality:_advance(tCurr, inFld, outFld)
   self.advanceFunc(tCurr, inFld, outFld)
end

function SpitzerCollisionality:_advanceOnDevice(tCurr, inFld, outFld)
   self.advanceOnDeviceFunc(tCurr, inFld, outFld)
end

return SpitzerCollisionality
