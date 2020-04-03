-- Gkyl ------------------------------------------------------------------------
--
-- Updater to (directly) solve generalized Poisson equation 
--
--      nabla_i D^ij nabla_j phi  = - rho
--
-- in 2D with RDG discretization.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

ffi.cdef[[
  typedef struct DiscontPoisson DiscontPoisson;
  DiscontPoisson* new_DiscontPoisson(int ncells[3], int ndim, int nbasis, int nnonzero, int polyOrder, bool writeMatrix);
  void delete_DiscontPoisson(DiscontPoisson* f);

  void discontPoisson_pushTriplet(DiscontPoisson* f, int idxK, int idxL, double val);
  void discontPoisson_constructStiffMatrix(DiscontPoisson* f);
  void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src, double* srcMod);
  void discontPoisson_solve(DiscontPoisson* f);
  void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol);
]]

-- DG Poisson solver updater object.
local DiscontGenPoisson = Proto(UpdaterBase)

function DiscontGenPoisson:init(tbl)
   DiscontGenPoisson.super.init(self, tbl)

   self.grid = assert(tbl.onGrid, "Updater.DiscontGenPoisson: Must provide grid object using 'onGrid'")
   self.basis = assert(tbl.basis, "Updater.DiscontGenPoisson: Must specify basis functions to use using 'basis'")

   assert(self.grid:ndim() == self.basis:ndim(),
          "Dimensions of basis and grid must match")
   
   self.ndim = self.grid:ndim()
   self.nbasis = self.basis:numBasis()
   local polyOrder = self.basis:polyOrder()
   
   self.ncell = ffi.new("int[3]")
   self.dx = Lin.Vec(3)    -- Limited to uniform grids for now.
   for d = 1,self.ndim do
      self.ncell[d-1] = self.grid:numCells(d)
      self.dx[d] = self.grid:dx(d)
   end

   local writeMatrix = xsys.pickBool(tbl.writeMatrix, false)

   -- Read the boundary conditions in.
   assert(#tbl.bcLower == self.ndim, "Updater.DiscontGenPoisson: Must provide lower boundary conditions for all the dimesions using 'bcLower'")
   assert(#tbl.bcUpper == self.ndim, "Updater.DiscontGenPoisson: Must provide upper boundary conditions for all the dimesions using 'bcUpper'")
   self.bcLower = tbl.bcLower
   self.bcUpper = tbl.bcUpper

   local basisNm = ''
   if self.ndim > 1 and polyOrder > 1 then
      if self.basis:id() == 'serendipity' then
         basisNm = 'Ser'
      elseif self.basis:id() == 'tensor' then
         basisNm = 'Tensor'
      end
   end
      
   self.Dxx = tbl.Dxx
   self.Dyy = tbl.Dyy
   self.Dxy = tbl.Dxy
   
   self._first = true
   self._matrixFn = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_T = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_T_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_B = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_B_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_L = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_L_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_R = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_R_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_TL = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_T_L_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_TR = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_T_R_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_BL = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_B_L_%dp", basisNm, self.ndim, polyOrder))
   self._matrixFn_BR = require(string.format("Updater.discontGenPoissonData.discontGenPoisson%sStencil%dD_B_R_%dp", basisNm, self.ndim, polyOrder))

   self.nnonzero = self.nbasis^2
   self.poisson = ffiC.new_DiscontPoisson(self.ncell, self.ndim, self.nbasis,
                                          self.nnonzero, polyOrder, writeMatrix)

   return self
end

function DiscontGenPoisson:buildStiffMatrix(phi)
   local ndim = self.ndim

   local localRange = phi:localRange()
   local indexer = phi:genIndexer()
   local lower, upper = {}, {}
   for d = 1,ndim do
      lower[d] = localRange:lower(d)
      upper[d] = localRange:upper(d)
   end
   lower[ndim+1] = 1
   upper[ndim+1] = self.nbasis
   local stiffMatrixRange = Range.Range(lower, upper)
   local stiffMatrixIndexer = Range.makeRowMajorGenIndexer(stiffMatrixRange)

   lower, upper = {}, {}
   for d = 1,ndim do
      lower[d] = -1
      upper[d] = 1
   end
   local stencilRange = Range.Range(lower, upper)
   local stencilIndexer = Range.makeRowMajorGenIndexer(stencilRange)

   local idxsExtK, idxsExtL = Lin.Vec(ndim+1), Lin.Vec(ndim+1)
   local cnt, val = 1, 0.0
   local idxK, idxL = 0, 0

   local SM = {}
   
   for idxs in localRange:colMajorIter() do
      local DxxPtr = self.Dxx:get(indexer(idxs))
      local DyyPtr = self.Dyy:get(indexer(idxs))
      local DxyPtr = self.Dxy:get(indexer(idxs))

      local SM = {}
      if idxs[1] == localRange:lower(1) and idxs[2] == localRange:lower(2) then
         SM = self._matrixFn_BL(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcLower[1].D,
                                self.bcLower[1].N,
                                self.bcLower[1].val,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[1] == localRange:lower(1) and idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_TL(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcLower[1].D,
                                self.bcLower[1].N,
                                self.bcLower[1].val,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:lower(2) then
         SM = self._matrixFn_BR(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_TR(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      elseif idxs[1] == localRange:lower(1) then
         SM = self._matrixFn_L(self.dx, DxxPtr, DyyPtr, DxyPtr,
                               self.bcLower[1].D,
                               self.bcLower[1].N,
                               self.bcLower[1].val,
                               0, 0, 0)
      elseif idxs[1] == localRange:upper(1) then
         SM = self._matrixFn_R(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                0, 0, 0)
      elseif idxs[2] == localRange:lower(2)then
         SM = self._matrixFn_B(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                0, 0, 0,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_T(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                0, 0, 0,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      else
         SM = self._matrixFn(self.dx, DxxPtr, DyyPtr, DxyPtr,
                             0, 0, 0, 0, 0, 0)
      end

      for d = 1,ndim do
         idxsExtK[d] = idxs[d]
      end
      for k = 1,self.nbasis do
         idxsExtK[ndim+1] = k
         idxK = stiffMatrixIndexer(idxsExtK)
         for l = 1,self.basis:numBasis() do
            for stencilIdx in stencilRange:colMajorIter() do
               for d = 1,ndim do
                  idxsExtL[d] = idxs[d] + stencilIdx[d]
               end
               idxsExtL[ndim+1] = l
               idxL = stiffMatrixIndexer(idxsExtL)
               val = SM[stencilIndexer(stencilIdx)][k][l]
               if math.abs(val) > 1e-14 then
                  ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
               end
            end        
         end
      end
   end
   ffiC.discontPoisson_constructStiffMatrix(self.poisson);
end

-- Advance method.
function DiscontGenPoisson:_advance(tCurr, inFld, outFld)
   local ndim = self.ndim

   local src = assert(inFld[1], "DiscontGenPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "DiscontGenPoisson.advance: Must specify an output field")

   local localRange = src:localRange()
   local lower, upper = {}, {}
   for d = 1,ndim do
      lower[d] = localRange:lower(d)
      upper[d] = localRange:upper(d)
   end
   lower[ndim+1] = 1
   upper[ndim+1] = self.nbasis
   local stiffMatrixRange = Range.Range(lower, upper)
   local stiffMatrixIndexer = Range.makeRowMajorGenIndexer(stiffMatrixRange)

   if self._first then
      self.srcIndexer = src:genIndexer()
      self.solIndexer = sol:genIndexer()

      -- Construct the stiffness matrix using Eigen.
      self:buildStiffMatrix(sol)
   end

   -- Pushing source to the Eigen matrix.
   local idxsExt = Lin.Vec(ndim+1)
   local srcMod = Lin.Vec(self.nbasis)
   for idxs in localRange:colMajorIter() do
      local DxxPtr = self.Dxx:get(self.srcIndexer(idxs))
      local DyyPtr = self.Dyy:get(self.srcIndexer(idxs))
      local DxyPtr = self.Dxy:get(self.srcIndexer(idxs))

      if idxs[1] == localRange:lower(1) and idxs[2] == localRange:lower(2) then
         SM = self._matrixFn_BL(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcLower[1].D,
                                self.bcLower[1].N,
                                self.bcLower[1].val,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[1] == localRange:lower(1) and idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_TL(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcLower[1].D,
                                self.bcLower[1].N,
                                self.bcLower[1].val,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:lower(2) then
         SM = self._matrixFn_BR(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_TR(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      elseif idxs[1] == localRange:lower(1) then
         SM = self._matrixFn_L(self.dx, DxxPtr, DyyPtr, DxyPtr,
                               self.bcLower[1].D,
                               self.bcLower[1].N,
                               self.bcLower[1].val,
                               0, 0, 0)
      elseif idxs[1] == localRange:upper(1) then
         SM = self._matrixFn_R(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                self.bcUpper[1].D,
                                self.bcUpper[1].N,
                                self.bcUpper[1].val,
                                0, 0, 0)
      elseif idxs[2] == localRange:lower(2)then
         SM = self._matrixFn_B(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                0, 0, 0,
                                self.bcLower[2].D,
                                self.bcLower[2].N,
                                self.bcLower[2].val)
      elseif idxs[2] == localRange:upper(2) then
         SM = self._matrixFn_T(self.dx, DxxPtr, DyyPtr, DxyPtr,
                                0, 0, 0,
                                self.bcUpper[2].D,
                                self.bcUpper[2].N,
                                self.bcUpper[2].val)
      else
         SM = self._matrixFn(self.dx, DxxPtr, DyyPtr, DxyPtr,
                             0, 0, 0, 0, 0, 0)
      end

      for k = 1,self.nbasis do
         srcMod[k] = -SM[10][k]
      end
      
      local srcPtr = src:get(self.srcIndexer(idxs))
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_pushSource(self.poisson, idx-1, srcPtr:data(), srcMod:data())
   end

   ffiC.discontPoisson_solve(self.poisson)

   for idxs in localRange:colMajorIter() do
      local solPtr = sol:get(self.solIndexer(idxs))
      local idxsExt = Lin.Vec(ndim+1)
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_getSolution(self.poisson, idx-1, solPtr:data())
   end

   self._first = false
end

function DiscontGenPoisson:delete()
   ffiC.delete_DiscontPoisson(self.poisson)
end

return DiscontGenPoisson
