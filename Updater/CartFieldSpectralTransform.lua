-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute transforms of DG data to the space spanned by
-- a global, spectral basis (e.g. Fourier, Hermite, Laguerre).
--
-- As a first step we will compute Hermite transforms of 1x1v data along the
-- velocity dimension. This entails:
-- 1) Projecting each of the 1D Hermite functions onto the DG basis.
-- 2) Allocating the Eigen objects. Will need one Nv*(p+1)xNv*(p+1) matrix
--    and its inverse, and one Nv*(p+1) vector for each of Nc surface basis elements.
-- 3) Initialize the matrix: Need to loop over velocity space and assign each
--    of the Hermite coefficients to the right matrix entry.
-- 4) Within a loop over configuration space:
--  4a. In a kernel, project f onto basis functions and assign coefficients of
--      the resulting 1D-like expansion to the corresponding Eigen vectors.
--  4b. Perform the Eigen multiplication A^{-1} b.
--  4c. Re-distribute data.
--
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local xsys           = require "xsys"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local DataStruct     = require "DataStruct"
local Basis          = require "Basis"
local Grid           = require "Grid"
local Hermite        = require "Lib.Hermite"
-- System libraries.
local ffi           = require "ffi"
local ffiC          = ffi.C

ffi.cdef[[
  typedef struct spectralTransform spectralTransform;
  spectralTransform* new_spectralTransform(const int nModes, const int nCells, const int pOrder, const int nSurfB);
  void assignLHSMatrixSer(spectralTransform *sTransObj, const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn);
  void getLHSMatrixInverse(spectralTransform *sTransObj);
  void assignRHSMatrix1x1vSer_P1OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P1OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P2OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P2OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void solveTransform(spectralTransform *sTransObj);
  void getSolution1x1vSer_P1(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);
  void getSolution1x1vSer_P2(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);
]]

-- Function to check if transform name is correct.
local function isTransformNameGood(nm)
--   if nm == "Fourier" or nm == "Hermite" or nm=="Laguerre" then
   if nm == "Hermite" then    -- Only support Hermite for now.
      return true
   end
   return false
end

-- Spetral transform updater object.
local CartFieldSpectralTransform = Proto(UpdaterBase)

function CartFieldSpectralTransform:init(tbl)
   CartFieldSpectralTransform.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.CartFieldSpectralTransform: Must provide grid object using 'onGrid'.")

   local weakBasis = assert(
      tbl.weakBasis, "Updater.CartFieldSpectralTransform: Must provide the weak basis object using 'weakBasis'.")
   
   -- If doing a kinetic simulation, also specify the configuration-space basis.
   -- Used to determine cDim and vDim, important for enforcing the assumption that
   -- Hermite and Laguerre transforms are only done in velocity space.
   local confBasis = tbl.confBasis

   -- Desired operation. The direction of the transform is given when calling the advance method.
   local op = assert(
      tbl.transforms, "Updater.CartFieldSpectralTransform: Must indicate desired transforms using the 'transforms' table.")

   -- Indicate in which direction to perform the transform.
   -- this is 1-indexed, and goes over all of (phase-space) directions.
   local opDir = assert(
      tbl.directions, "Updater.CartFieldSpectralTransform: Must indicate desired transform directions using the table 'directions'.")

   -- Dimension, name and polynomial order of weak basis, and number of basis functions.
   self._wDim      = weakBasis:ndim()
   self._basisID   = weakBasis:id()
   self._polyOrder = weakBasis:polyOrder()
   self._numBasis  = weakBasis:numBasis()
   if confBasis then
      self._cDim = confBasis:ndim()
      self._vDim = self._wDim - self._cDim
   end

   -- Number of spectral modes represented in each 'opDir' direction.
   self._numSpectralModes = {}
   for iDir = 1,#opDir do
      self._numSpectralModes[iDir] = 44 --(self._polyOrder+1)*self._onGrid:numCells(opDir[iDir])
   end

   -- Function to create basis functions.
   local function createBasis(nm, ndim, polyOrder)
      if nm == "serendipity" then
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      elseif nm == "maximal-order" then
         return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
      elseif nm == "tensor" then
         return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
      end
   end
   -- Number of surface basis elements. Equivalent to number of elements of one lower dimension.
   local surfBasis          = createBasis(self._basisID, self._wDim-1, self._polyOrder)
   self._numSurfBasis = surfBasis:numBasis()

   for iOp = 1,#op do
      if isTransformNameGood(op[iOp]) then
         if op[iOp]=="Hermite" or op[iOp]=="Laguerre" then
            assert(confBasis,
               "Updater.CartFieldSpectralTransform: Provide configuration-space basis with 'confBasis'")
            assert(opDir[iOp] > self._cDim,
               "Updater.CartFieldSpectralTransform: Hermite/Laguerre transforms only in velocity space." .. 
                " opDir must be greater than configuration space dimensionality.")
         end
         -- Select the right kernel to assign mass matrix elements, assign RHS matrix elements
         -- re-distribute solution from the resulting matrix to a Gkeyll CartField.
--         self._assignMatrixEs = TransformDecl.selectAssign(op, self._basisID, self._wDim, self._polyOrder)
      else
         assert(false, string.format(
           	"CartFieldSpectralTransform: transform must be one of Fourier, Hermite, Laguerre. Requested %s instead.", op[iOp]))
      end
   end

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   -- Create struct containing allocated SpectralTransform arrays for solving the linear global problem.
   self._transformObj = ffiC.new_spectralTransform(self._numSpectralModes[1], self._onGrid:numCells(opDir[1]),
                                                   self._polyOrder, self._numSurfBasis) 

   -- Create projection updater to project the spectral basis functions onto our basis. 
   -- Need to create a 1D grid to project the spectral basis onto.
   local basis1D       = createBasis(self._basisID, 1, self._polyOrder)
   local opDirPeriodic = {}
   if self._onGrid:isDirPeriodic(opDir[1]) then    -- Only one transform for now.
      opDirPeriodic[1] = 1
   end
   self._grid1D = Grid.RectCart {
      lower         = {self._onGrid:lower(opDir[1])},
      upper         = {self._onGrid:upper(opDir[1])},
      cells         = {self._onGrid:numCells(opDir[1])},
      periodicDirs  = {opDirPeriodic},
      -- These should match those of the parent grid.
--      decomposition = self.decomp,
--      mappings      = coordinateMap,
   }
   self._project = ProjectOnBasis {
      onGrid   = self._grid1D,
      basis    = basis1D,
      evaluate = function (t,xn)  -- The function is set below with the setFuncMethod.
                    return 1.0
                 end,
      projectOnGhosts = false,
   }
   -- Create a weak field for the basis element.
   local ghostCells = {0,0}
   if self.onGhosts then ghostCells = {1,1} end
   self._spectralBasisE = DataStruct.Field {
      onGrid        = self._grid1D,
      numComponents = basis1D:numBasis(),
      ghost         = ghostCells,
   }

--   self._spectralBasisEIndexer = self._spectralBasisE:genIndexer()
--   self._spectralBasisEPtr     = self._spectralBasisE:get(1)
--
--   self._spectralBasisERange   = self._spectralBasisE:localRange()
--   if self.onGhosts then    -- Extend range to config-space ghosts.
--      spectralBasisERange = spectralBasisERange:extendDir(1, spectralBasisE:lowerGhost(), spectralBasisE:upperGhost())
--   end
   self._isFirst   = true
   self._perpRange = {}    -- Ranges in directions perpendicular to opDir.

   -- Cell index and length.
   self.idx = Lin.IntVec(self._wDim)
   self.dx  = Lin.Vec(self._wDim)

   local iDir = 1      -- Only one transform for now.
   if op[iDir] == "Hermite" then

--      local spectralRangeDecomp = LinearDecomp.LinearDecompRange {
--         range = spectralBasisERange:selectFirst(1), numSplit = self._grid1D:numSharedProcs() }
--      local tId = self._grid1D:subGridSharedId()    -- Local thread ID.
--
--      for mH = 0, self._numSpectralModes[iDir]-1 do
--         -- For each Hermite basis function need to:
--         --   1) Project the Hermite basis function onto DG basis.
--         --   2) Assign its coefficients to the corresponding matrix element.
--
--         -- Project m-th hermite basis onto DG basis.
--         project:setFunc(function (t,xn) return Hermite.H(mH,xn[1]) end )
--         project:advance(0.0, {}, {spectralBasisE})
--
--         -- Loop over 1D space assigning elements to matrix for 1D spectral transform.
--         for sIdx in spectralRangeDecomp:rowMajorIter(tId) do
--
--            self._grid1D:setIndex(sIdx)
--            spectralBasisE:fill(spectralBasisEIndexer(sIdx), spectralBasisEPtr)
--
--            ffiC.assignLHSMatrixSer(self._transformObj, self._polyOrder, sIdx[iDir], mH, spectralBasisEPtr:data())
--
--         end
--      end
--
--      -- Invert the left-side mass matrix.
--      ffiC.getLHSMatrixInverse(self._transformObj)

-- Only Hermite is supported for now.
--   elseif op[iDir] == "Laguerre" then
--   elseif op[iDir] == "Fourier" then
   end

end

-- Advance method.
function CartFieldSpectralTransform:_advance(tCurr, inFld, outFld)
   local grid             = self._onGrid
   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local dgFld       = inFld[1]
   local spectralFld = outFld[1]

   local phaseRange = dgFld:localRange()
--   if self.onGhosts then    -- Extend range to config-space ghosts.
--      for dir = 1, cDim do
--         phaseRange = phaseRange:extendDir(dir, dgFld:lowerGhost(), dgFld:upperGhost())
--      end
--   end

   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = grid:subGridSharedId()    -- Local thread ID.

   local dgFldIndexer       = dgFld:genIndexer()
   local spectralFldIndexer = spectralFld:genIndexer()

   local dgFldPtr       = dgFld:get(1)
   local spectralFldPtr = spectralFld:get(1)

   local spectralBasisEIndexer = self._spectralBasisE:genIndexer()
   local spectralBasisEPtr     = self._spectralBasisE:get(1)

   local spectralBasisERange   = self._spectralBasisE:localRange()

   local spectralRangeDecomp = LinearDecomp.LinearDecompRange {
      range = spectralBasisERange:selectFirst(1), numSplit = self._grid1D:numSharedProcs() }
   local tId = self._grid1D:subGridSharedId()    -- Local thread ID.

   -- Outer loop is threaded and over configuration space.
   -- This will later be over perpendicular dimensions.
   for cIdx in confRangeDecomp:rowMajorIter(tId) do

      for mH = 0, self._numSpectralModes[1]-1 do
         -- For each Hermite basis function need to:
         --   1) Project the Hermite basis function onto DG basis.
         --   2) Assign its coefficients to the corresponding matrix element.

         -- Project m-th hermite basis onto DG basis.
         self._project:setFunc(function (t,xn) return Hermite.H(mH,xn[1]) end )
         self._project:advance(0.0, {}, {self._spectralBasisE})

         -- Loop over 1D space assigning elements to matrix for 1D spectral transform.
         for sIdx in spectralRangeDecomp:rowMajorIter(tId) do

            self._grid1D:setIndex(sIdx)
            self._spectralBasisE:fill(spectralBasisEIndexer(sIdx), spectralBasisEPtr)

            ffiC.assignLHSMatrixSer(self._transformObj, self._polyOrder, sIdx[1], mH, spectralBasisEPtr:data())

         end
      end

      -- Invert the left-side mass matrix.
      ffiC.getLHSMatrixInverse(self._transformObj)

      cIdx:copyInto(self.idx)

      for vDir = 1, vDim do

         -- Lower/upper bounds in direction 'vDir': cell indices.
         local dirLoIdx, dirUpIdx = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)

         if self._isFirst then
            -- Restricted velocity range.
            -- Velocity integral in m0Star does not include last cell.
            self._perpRange[vDir] = phaseRange
            for cd = 1, cDim do
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
            end
            self._perpRange[vDir] = self._perpRange[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
         end
         local perpRange = self._perpRange[vDir]

         for vPerpIdx in perpRange:rowMajorIter() do
            vPerpIdx:copyInto(self.idx)
            for d = 1, cDim do self.idx[d] = cIdx[d] end

            -- Loop over transform dimension and assign the right-side (source) matrix.
            for i = dirLoIdx, dirUpIdx do     -- This loop is over cells.
               self.idx[cDim+vDir] = i

               grid:setIndex(self.idx)

               dgFld:fill(dgFldIndexer(self.idx), dgFldPtr)
               ffiC.assignRHSMatrix1x1vSer_P1OpDir2(self._transformObj, i, dgFldPtr:data())
            end

            -- Solve the linear transform problem.
            ffiC.solveTransform(self._transformObj)

            -- Loop over transform dimension and re-distribute solution into a CartField.
            --for i = dirLoIdx, dirUpIdx do     -- This loop is over cells.
            for i = 1, self._numSpectralModes[1]/(self._numBasis-self._numSurfBasis) do     -- This loop is over cells.
               self.idx[cDim+vDir] = i

               grid:setIndex(self.idx)

               spectralFld:fill(spectralFldIndexer(self.idx), spectralFldPtr)
               ffiC.getSolution1x1vSer_P1(self._transformObj, i, spectralFldPtr:data())
            end    -- end loop over transform dimension.
         end    -- end loop over perpendicular dimensions.
      end    -- end loop over velocity dimensions.
   end    -- end configuration space loop.
   self._isFirst = false
end

return CartFieldSpectralTransform
