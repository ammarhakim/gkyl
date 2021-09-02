-- Gkyl ------------------------------------------------------------------------
--
-- Apply twist shift operation in y to a Cartesian field (e.g. gyrokinetics),
-- using weak equality-based interpolation. Essentially we are computing the
-- left side of g(x,y) = f(x,y-S(x)).
--
-- Current limitations:
--    1) The shift has to be monotonic.
--    2) The shift can't be constant and a multiple of y-cell length anywhere.
--    3) Apply twist-shift to one 2D field and output to another 2D field, or
--       apply the twist-shift to a 3D field (filling the ghost cells).
--
-- Notes:
--    a] Need to figure out how to do this more accurately. In 2D p=1 the last
--       DG coefficient still seems wrong. I now suspect that this is because I
--       have not been projecting in the space of the shifted region, but rather
--       in the space of the donor cell that overlaps with the shifted region.
--    b] Need to experiment more and check robustness.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local Grid         = require "Grid"
local Basis        = require "Basis"
local DataStruct   = require "DataStruct"
local EvOnNodesUpd = require "Updater.EvalOnNodes"
local Range        = require "Lib.Range"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local math         = require "sci.math"  -- For sign function.

local TwistShiftDecl = require "Updater.twistShiftData.TwistShiftModDecl"
local tsFun          = require "Updater.twistShiftData.TwistShiftFun"

local TwistShiftBC = Proto(UpdaterBase)

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

function TwistShiftBC:init(tbl)
   TwistShiftBC.super.init(self, tbl) -- Setup base object.

   self.grid = assert(
      tbl.onGrid, "Updater.TwistShift: Must provide grid to interpolate to using 'onGrid'.")
   self.basis = assert(
      tbl.basis, "Updater.TwistShift: Must provide the basis of the fields using 'basis'.")

   self.confBasis = tbl.confBasis
   if self.confBasis then
      self.cDim = self.confBasis:ndim()
      self.vDim = self.basis:ndim() - self.cDim
   else
      self.cDim, self.vDim = self.basis:ndim(), 0
   end

   if self.cDim == 3 then
      self.zEdge = assert(
         tbl.edge, "Updater.TwistShift: Must indicate which z-edge to compute (lower/upper) using 'edge'.")
      assert(self.zEdge=="lower" or self.zEdge=="upper", "Updater.TwistShift: 'edge' must be lower or upper.")
   end

   local yShFunc = assert(
      tbl.yShiftFunc, "Updater.TwistShift: Must provide the y-shift function using 'yShiftFunc'.")
   local yShPolyOrder = assert(
      tbl.yShiftPolyOrder, "Updater.TwistShift: Must provide the polyOrder for projecting y-shift using 'yShiftPolyOrder'.")

   -- Project the y-shift function (of x) onto a 1D grid/basis.
   local yShGridIngr = self.grid:childGrid({1})
   local yShGrid = Grid.RectCart {
      lower = yShGridIngr.lower,  periodicDirs  = yShGridIngr.periodicDirs,
      upper = yShGridIngr.upper,  decomposition = yShGridIngr.decomposition,
      cells = yShGridIngr.cells,
   }
   local yShBasis = createBasis(yShGrid:ndim(), yShPolyOrder, self.basis:id())
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
   --print(Mpi.Comm_rank(grid:commSet().sharedComm), "before TS init")

   -- Call function computing the donor cells.
   self.doCells = tsFun.getDonors(self.grid, self.yShFld, yShBasis)

   -- Allocate matrices that multiply each donor cell to compute its contribution to a target cell.
   -- Also allocate a temp vector and matrix used in the mat-vec multiply.
   self.matVec = tsFun.matVec_alloc(self.yShFld, self.doCells, self.basis)

   -- Select kernels that assign matrices and later, in the :_advance method, do mat-vec multiplies.
   tsFun.selectTwistShiftKernels(self.cDim, self.vDim, self.basis:id(), self.basis:polyOrder(), yShPolyOrder)

   -- Pre-compute matrices using weak equalities between donor and target fields.
   tsFun.preCalcMat(self.grid, self.yShFld, self.doCells, self.matVec)

   self.tsMatVecMult = TwistShiftDecl.selectTwistShiftMatVecMult(self.cDim, self.vDim, self.basis:id(), self.basis:polyOrder())

   self.idxDoP = Lin.IntVec(self.cDim+self.vDim)

   self.isFirst = true   -- Will be reset the first time _advance() is called.

end

function TwistShiftBC:_advance(tCurr, inFld, outFld)
   local fldDo, fldTar = inFld[1], outFld[1]

   local cDim, vDim = self.cDim, self.vDim

   local oneFld = false
   if fldDo==nil then
      assert(cDim==3, "Updater.TwistShift: performing twist-shift to a single field is only available for 3D fields.")
      oneFld = true
      fldDo  = outFld[1]
   else
      assert(cDim==2, "Updater.TwistShift: donor and target fields are separate, so they must be 2D fields.")
   end

   local indexer = fldTar:genIndexer()
   local fldDoItr, fldTarItr = fldDo:get(1), fldTar:get(1)

   if oneFld then
 
      -- fetch skin cells from opposite edge in z (donors) 
      -- and copy to ghost cells (targets)
      local fldGrid = fldTar:grid()
      local periodicDirs = fldGrid:getPeriodicDirs()
      local modifiedPeriodicDirs = {3}
      fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
      fldTar:sync()
      fldGrid:setPeriodicDirs(periodicDirs)

      local global = fldTar:globalRange()

      if self.isFirst then
         local getGhostRange = function(globalIn, globalExtIn, dir, edge)
            local lv, uv = globalIn:lowerAsVec(), globalIn:upperAsVec()
            if edge == "lower" then
               lv[dir] = globalIn:lower(dir)-1
               uv[dir] = lv[dir]
            else
               lv[dir] = globalIn:upper(dir)+1
               uv[dir] = lv[dir]
            end
            return Range.Range(lv, uv)
         end
         local globalExt     = fldTar:globalExtRange()
         local localExtRange = fldTar:localExtRange()
         self.ghostRng = localExtRange:intersect(getGhostRange(global, globalExt, 3, self.zEdge))
         -- Decompose ghost region into threads.
         self.ghostRangeDecomp = LinearDecomp.LinearDecompRange {
            range = self.ghostRng, numSplit = self.grid:numSharedProcs() }

         self.isFirst = false
      end

      local tId = self.grid:subGridSharedId() -- Local thread ID.

      for idxTar in self.ghostRangeDecomp:rowMajorIter(tId) do

         fldTar:fill(indexer(idxTar), fldTarItr)
         -- Zero out target cell before operation.
         for i = 1, self.basis:numBasis() do fldTarItr[i] = 0. end

         local doCellsC = self.doCells[idxTar[1]][idxTar[2]]

         idxTar:copyInto(self.idxDoP)
         -- get z index of skin cells (donors) that will be shift-copied to ghost cells (targets)
         self.idxDoP[3] = self.zEdge=="lower" and global:upper(3) or global:lower(3)
--         print("idxTar = ",idxTar[1],idxTar[2],idxTar[3],idxTar[4],idxTar[5])

         for mI = 1, #doCellsC do

            local idxDo2D = doCellsC[mI]
            self.idxDoP[1], self.idxDoP[2] = idxDo2D[1], idxDo2D[2] 
--            print("   from idxDo = ",self.idxDoP[1],self.idxDoP[2],self.idxDoP[3],self.idxDoP[4],self.idxDoP[5])

            fldDo:fill(indexer(self.idxDoP), fldDoItr)

            -- Matrix-vec multiply to compute the contribution of each donor cell to a target cell..
            self.tsMatVecMult(self.matVec, idxTar[1], mI, fldDoItr:data(), fldTarItr:data())
         end
      end
      
   else   -- Only for 2D fields.
      local localRange = fldTar:localRange()

      for idxTar in localRange:rowMajorIter() do

         fldTar:fill(indexer(idxTar), fldTarItr)
         -- Zero out target cell before operation.
         for i = 1, self.basis:numBasis() do fldTarItr[i] = 0. end

         local doCellsC = self.doCells[idxTar[1]][idxTar[2]]

         idxTar:copyInto(self.idxDoP)
--         print("idxTar = ",idxTar[1],idxTar[2])

         for mI = 1, #doCellsC do

            local idxDo2D = doCellsC[mI]
--            if idxTar[2]==10 and idxDo2D[2]==6 then
            self.idxDoP[1], self.idxDoP[2] = idxDo2D[1], idxDo2D[2] 
--            print("   from idxDo = ",self.idxDoP[1],self.idxDoP[2])

            fldDo:fill(indexer(self.idxDoP), fldDoItr)

            -- Matrix-vec multiply to compute the contribution of each donor cell to a target cell..
            self.tsMatVecMult(self.matVec, idxTar[1], mI, fldDoItr:data(), fldTarItr:data())
--            end
         end
      end
   end

end

return TwistShiftBC
