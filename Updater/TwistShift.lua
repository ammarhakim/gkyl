-- Gkyl ------------------------------------------------------------------------
--
-- Apply twist shift operation in y to a Cartesian field (e.g. gyrokinetics),
-- using weak equality-based interpolation.
--
-- Current limitations:
--    1) The shift has to be monotonic.
--    2) The shift can't be constant and a multiple of y-cell length anywhere.
--    3) Apply twist-shift to one 2D field and output to another 2D field, or
--       apply the twist-shift to a 3D field (filling the ghost cells).
--
-- Notes:
--    a] Need to figure out how to do this more accurately. In 2D p=1 the last
--       DG coefficient still seems wrong.
--    b] Need to experiment more and check robustness.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local Grid           = require "Grid"
local Basis          = require "Basis"
local DataStruct     = require "DataStruct"
local EvOnNodesUpd   = require "Updater.EvalOnNodes"

local TwistShiftDecl = require "Updater.twistShiftData.TwistShiftModDecl"
local tsFun          = require "Updater.twistShiftData.TwistShiftFun"

local TwistShift = Proto(UpdaterBase)

local function createBasis(dim, pOrder, bKind)
   local basis
   if (bKind=="serendipity") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="maximal-order") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

function TwistShift:init(tbl)
   TwistShift.super.init(self, tbl) -- Setup base object.

   self.grid = assert(
      tbl.onGrid, "Updater.TwistShift: Must provide grid to interpolate to using 'onGrid'.")
   local basis = assert(
      tbl.basis, "Updater.TwistShift: Must provide the basis of the fields using 'basis'.")

   local confBasis = tbl.confBasis
   local cDim, vDim
   if confBasis then
      cDim = confBasis:ndim()
      vDim = basis:ndim() - cDim
   else
      cDim = basis:ndim()
      vDim = 0
   end

   local yShFunc = assert(
      tbl.yShiftFunc, "Updater.TwistShift: Must provide the y-shift function using 'yShiftFunc'.")
   local yShPolyOrder = assert(
      tbl.yShiftPolyOrder, "Updater.TwistShift: Must provide the polyOrder for projecting y-shift using 'yShiftPolyOrder'.")

   -- Project the y-shift function onto a 1D grid/basis.
   local yShGridIngr = self.grid:childGrid({1})
   local yShGrid = Grid.RectCart {
      lower = yShGridIngr.lower,
      upper = yShGridIngr.upper,
      cells = yShGridIngr.cells,
      periodicDirs  = yShGridIngr.periodicDirs,
      decomposition = yShGridIngr.decomposition,
   }
   local yShBasis = createBasis(yShGrid:ndim(), yShPolyOrder, basis:id())
   self.yShFld = DataStruct.Field {
      onGrid        = yShGrid,
      numComponents = yShBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = {polyOrder = yShBasis:polyOrder(), basisType = yShBasis:id()},
   }
   local projUpd = EvOnNodesUpd {
      onGrid   = yShGrid,
      basis    = yShBasis,
      evaluate = yShFunc,
   }
   projUpd:advance(0., {}, {self.yShFld})

   -- Set some constant variables needed across tsFun functions.
   tsFun.set_yShiftF(yShFunc)
   tsFun.set_domLim(self.grid)
   tsFun.set_dx(self.grid)
   tsFun.set_projData(yShPolyOrder)

   -- Call function computing the donor cells.
   self.doCells = tsFun.getDonors(self.grid, self.yShFld, yShBasis)

   -- Allocate matrices that multiply each donor cell to compute its contribution to a target cell.
   -- Also allocate a temp vector and matrix used in the mat-vec multiply.
   self.matVec = tsFun.matVec_alloc(self.yShFld, self.doCells, basis)

   -- Select kernels that assign matrices and later, in the :_advance method, do mat-vec multiplies.
   tsFun.selectTwistShiftKernels(cDim, vDim, basis:id(), basis:polyOrder(), yShPolyOrder)

   -- Pre-compute matrices using weak equalities between donor and target fields.
   tsFun.preCalcMat(self.grid, self.yShFld, self.doCells, self.matVec)

   self.tsMatVecMult = TwistShiftDecl.selectTwistShiftMatVecMult(cDim, vDim, basis:id(), basis:polyOrder())

end

function TwistShift:_advance(tCurr, inFld, outFld)
   local fldDo, fldTar = inFld[1], outFld[1]

   local localRange = fldTar:localRange()

   local indexer             = fldTar:genIndexer()
   local fldDoItr, fldTarItr = fldDo:get(1), fldTar:get(1)

   for idxTar in localRange:rowMajorIter() do

      fldTar:fill(indexer(idxTar), fldTarItr)

      local doCellsC = self.doCells[idxTar[1]][idxTar[2]]

      for mI = 1, #doCellsC do

         fldDo:fill(indexer(doCellsC[mI]), fldDoItr)

         -- Matrix-vec multiply to compute the contribution of each donor cell to a target cell..
         self.tsMatVecMult(self.matVec, idxTar[1], mI, fldDoItr:data(), fldTarItr:data())
      end
   end

end

return TwistShift
