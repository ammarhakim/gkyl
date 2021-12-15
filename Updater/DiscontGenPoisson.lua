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

-- Gkyl libraries
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

ffi.cdef[[
  typedef struct DiscontPoisson DiscontPoisson;
  DiscontPoisson* new_DiscontPoisson(const char* outPrefix,
     int ncells[3], int ndim, int nbasis, int nnonzero, int polyOrder, bool writeMatrix);
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
   assert(self.grid:ndim() == 2, "DiscontGenPoisson is currently implemented only for 2D")
   
   self.ndim = self.grid:ndim()
   self.numBasis = self.basis:numBasis()
   local polyOrder = self.basis:polyOrder()

   self.ncell = ffi.new("int[3]")
   self.dx = Lin.Vec(3) -- uniform grids for now
   for d = 1,self.ndim do
      self.ncell[d-1] = self.grid:numCells(d)
      self.dx[d] = self.grid:dx(d)
   end

   local writeMatrix = xsys.pickBool(tbl.writeMatrix, false)

   -- Read boundary conditions
   assert(#tbl.bcLower == self.ndim, "Updater.DiscontGenPoisson: Must provide lower boundary conditions for all the dimesions using 'bcLower'")
   assert(#tbl.bcUpper == self.ndim, "Updater.DiscontGenPoisson: Must provide upper boundary conditions for all the dimesions using 'bcUpper'")
   self.bcLower, self.bcUpper = {}, {}
   for d = 1, self.ndim do
      if tbl.bcLower[d].T == "D" then
         self.bcLower[d] = {D=1, N=0, val=tbl.bcLower[d].V}
      elseif tbl.bcLower[d].T == "N" then
         self.bcLower[d] = {D=0, N=1, val=tbl.bcLower[d].V}
      else
         self.bcLower[d] = tbl.bcLower[d]
      end
      if tbl.bcUpper[d].T == "D" then
         self.bcUpper[d] = {D=1, N=0, val=tbl.bcUpper[d].V}
      elseif tbl.bcUpper[d].T == "N" then
         self.bcUpper[d] = {D=0, N=1, val=tbl.bcUpper[d].V}
      else
         self.bcUpper[d] = tbl.bcUpper[d]
      end
   end

   -- -- Project on basis functions
   -- local basis1D = Basis.CartModalSerendipity {
   --    ndim = 1,
   --    polyOrder = polyOrder
   -- }
   -- self.numBasis1D = basis1D:numBasis()
   -- self.numQuad = tbl.numQuad and tbl.numQuad or self.basis:polyOrder()+1 -- Number of quadrature points in each direction
   -- self.ordinates = GaussQuadRules.ordinates[self.numQuad]
   -- self.weights = GaussQuadRules.weights[self.numQuad]
   -- self.basisAtOrdinates = Lin.Mat(self.numQuad, self.numBasis1D)
   -- for n = 1, self.numQuad do
   --    basis1D:evalBasis({self.ordinates[n]}, self.basisAtOrdinates[n])
   -- end


   -- Read the diffusion tensor function
   self.Dxx = tbl.Dxx
   self.Dyy = tbl.Dyy
   self.Dxy = tbl.Dxy
   local function getField()
      return DataStruct.Field {
         onGrid = self.grid,
         numComponents = self.basis:numBasis(),
         ghost = {1, 1},
         metaData = {
            polyOrder = self.basis:polyOrder(),
            basisType = self.basis:id(),
         },
      }
   end

   -- local projectUpd = ProjectOnBasis {
   --    onGrid = self.grid,
   --    basis = self.basis,
   --    numQuad = self.numQuad,
   --    evaluate = function(t,xn) return 0. end,
   --    onGhosts = true,
   -- }

   -- local evaluateUpd = EvaluateOnBasis {
   --    onGrid = self.grid,
   --    basis = self.basis,
   --    numQuad = self.numQuad,
   --    evaluate = function(t,xn) return 0. end,
   --    onGhosts = true,
   -- }
   
   -- self.Dxx = getField()
   -- projectUpd:setFunc(function(t,xn) return self.DxxFn(t,xn) end)
   -- projectUpd:advance(0.0, {}, {self.Dxx})

   -- self.Dxy = getField()
   -- projectUpd:setFunc(function(t,xn) return self.DxyFn(t,xn) end)
   -- projectUpd:advance(0.0, {}, {self.Dxy})

   -- self.Dyy = getField()
   -- projectUpd:setFunc(function(t,xn) return self.DyyFn(t,xn) end)
   -- projectUpd:advance(0.0, {}, {self.Dyy})

   -- Load the kernels
   if self.basis:id() == 'serendipity' then
      basisNm = 'Ser'
   elseif self.basis:id() == 'tensor' then
      basisNm = 'Tensor'
   end

   -- create indexer for use in indexing stiffness matrix
   local lower, upper = {}, {}
   for d = 1, self.ndim do
      lower[d] = 1
      upper[d] = self.grid:numCells(d)
   end
   lower[self.ndim+1] = 1
   upper[self.ndim+1] = self.numBasis
   local stiffMatrixRange = Range.Range(lower, upper)
   self.stiffMatrixIndexer = Range.makeRowMajorGenIndexer(stiffMatrixRange)

   self._matrixFn = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_T = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_T_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_B = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_B_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_L = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_L_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_R = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_R_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_TL = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_T_L_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_TR = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_T_R_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_BL = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_B_L_p%d",
                    basisNm:lower(), self.ndim, polyOrder))
   self._matrixFn_BR = require(
      string.format("Updater.discontGenPoissonData.dg_poisson_anisotropic_%s_%dx_B_R_p%d",
                    basisNm:lower(), self.ndim, polyOrder))

   self.nnonzero = self.numBasis^2
   self.poisson = ffiC.new_DiscontPoisson(GKYL_OUT_PREFIX, self.ncell, self.ndim, self.numBasis,
                                          self.nnonzero, polyOrder, writeMatrix)

   self:buildStiffMatrix()

   return self
end

-- Get list of blocks in stiffness matrix for given index
-- location. List is returned as a set of 10 blocks:
--
-- { BL, L, TL, B, C, T, BR, R, TR, BOUNDARY }
--
-- BOUNDARY are boundary term modifications to RHS due to boundary
-- conditions
function DiscontGenPoisson:getBlock(idxs)
   local indexer = self.Dxx:genIndexer()
   local localRange = self.Dxx:localRange()

   local idxsT = {}
   local idxsB = {}
   local idxsL = {}
   local idxsR = {}

   idxsT[1], idxsT[2] = idxs[1], idxs[2]+1
   idxsB[1], idxsB[2] = idxs[1], idxs[2]-1
   idxsL[1], idxsL[2] = idxs[1]-1, idxs[2]
   idxsR[1], idxsR[2] = idxs[1]+1, idxs[2]
   
   -- as we only use two-cell recovery for Dij only values in five
   -- face neighbors of cell idxs are needed
   local DxxCPtr = self.Dxx:get(indexer(idxs))
   local DyyCPtr = self.Dyy:get(indexer(idxs))
   local DxyCPtr = self.Dxy:get(indexer(idxs))

   local DxxTPtr = self.Dxx:get(indexer(idxsT))
   local DyyTPtr = self.Dyy:get(indexer(idxsT))
   local DxyTPtr = self.Dxy:get(indexer(idxsT))

   local DxxBPtr = self.Dxx:get(indexer(idxsB))
   local DyyBPtr = self.Dyy:get(indexer(idxsB))
   local DxyBPtr = self.Dxy:get(indexer(idxsB))

   local DxxLPtr = self.Dxx:get(indexer(idxsL))
   local DyyLPtr = self.Dyy:get(indexer(idxsL))
   local DxyLPtr = self.Dxy:get(indexer(idxsL))

   local DxxRPtr = self.Dxx:get(indexer(idxsR))
   local DyyRPtr = self.Dyy:get(indexer(idxsR))
   local DxyRPtr = self.Dxy:get(indexer(idxsR))

   -- local numBasis = self.numBasis1D
   -- local DxxLF, DxyLF = Lin.Vec(numBasis), Lin.Vec(numBasis)
   -- local DxxRF, DxyRF = Lin.Vec(numBasis), Lin.Vec(numBasis)
   -- local DyyBF, DxyBF = Lin.Vec(numBasis), Lin.Vec(numBasis)
   -- local DyyTF, DxyTF = Lin.Vec(numBasis), Lin.Vec(numBasis)

   local xc = Lin.Vec(self.ndim)
   local dx = Lin.Vec(self.ndim)
   self.grid:setIndex(idxs)
   self.grid:getDx(dx)
   self.grid:cellCenter(xc)
   -- local z = Lin.Vec(2)
   -- for k = 1,numBasis do
   --    -- Left face
   --    DxxLF[k], DxyLF[k] = 0.0, 0.0
   --    z[1] = xc[1] - dx[1]/2
   --    for i = 1,self.numQuad do
   --       z[2] = xc[2] + self.ordinates[i]*dx[2]/2
   --       DxxLF[k] = DxxLF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxxFn(0, z)
   --       DxyLF[k] = DxyLF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxyFn(0, z)
   --    end

   --    -- Right face
   --    DxxRF[k], DxyRF[k] = 0.0, 0.0
   --    z[1] = xc[1] + dx[1]/2
   --    for i = 1,self.numQuad do
   --       z[2] = xc[2] + self.ordinates[i]*dx[2]/2
   --       DxxRF[k] = DxxRF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxxFn(0, z)
   --       DxyRF[k] = DxyRF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxyFn(0, z)
   --    end

   --    -- Bottom face
   --    DyyBF[k], DxyBF[k] = 0.0, 0.0
   --    z[2] = xc[2] - dx[2]/2
   --    for i = 1,self.numQuad do
   --       z[1] = xc[1] + self.ordinates[i]*dx[1]/2
   --       DyyBF[k] = DyyBF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DyyFn(0, z)
   --       DxyBF[k] = DxyBF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxyFn(0, z)
   --    end

   --    -- Top face
   --    DyyTF[k], DxyTF[k] = 0.0, 0.0
   --    z[2] = xc[2] + dx[2]/2
   --    for i = 1,self.numQuad do
   --       z[1] = xc[1] + self.ordinates[i]*dx[1]/2
   --       DyyTF[k] = DyyTF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DyyFn(0, z)
   --       DxyTF[k] = DxyTF[k] + self.weights[i]*self.basisAtOrdinates[i][k]*self.DxyFn(0, z)
   --    end
   -- end

   -- Fetch blocks base on cell index. We need special blocks for each
   -- corner cell, skin cells and interior cells
   local SM = nil
   if idxs[1] == localRange:lower(1) and idxs[2] == localRange:lower(2) then
      SM = self._matrixFn_BL(self.dx,
                             DxxCPtr, DyyCPtr, DxyCPtr,
                             DxxLPtr, DyyLPtr, DxyLPtr,
                             DxxRPtr, DyyRPtr, DxyRPtr,
                             DxxBPtr, DyyBPtr, DxyBPtr,
                             DxxTPtr, DyyTPtr, DxyTPtr,
                             self.bcLower[1].D,
                             self.bcLower[1].N,
                             self.bcLower[1].val,
                             self.bcLower[2].D,
                             self.bcLower[2].N,
                             self.bcLower[2].val)
   elseif idxs[1] == localRange:lower(1) and idxs[2] == localRange:upper(2) then
      SM = self._matrixFn_TL(self.dx,
                             DxxCPtr, DyyCPtr, DxyCPtr,
                             DxxLPtr, DyyLPtr, DxyLPtr,
                             DxxRPtr, DyyRPtr, DxyRPtr,
                             DxxBPtr, DyyBPtr, DxyBPtr,
                             DxxTPtr, DyyTPtr, DxyTPtr,
                             self.bcLower[1].D,
                             self.bcLower[1].N,
                             self.bcLower[1].val,
                             self.bcUpper[2].D,
                             self.bcUpper[2].N,
                             self.bcUpper[2].val)
   elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:lower(2) then
      SM = self._matrixFn_BR(self.dx,
                             DxxCPtr, DyyCPtr, DxyCPtr,
                             DxxLPtr, DyyLPtr, DxyLPtr,
                             DxxRPtr, DyyRPtr, DxyRPtr,
                             DxxBPtr, DyyBPtr, DxyBPtr,
                             DxxTPtr, DyyTPtr, DxyTPtr,
                             self.bcUpper[1].D,
                             self.bcUpper[1].N,
                             self.bcUpper[1].val,
                             self.bcLower[2].D,
                             self.bcLower[2].N,
                             self.bcLower[2].val)
   elseif idxs[1] == localRange:upper(1) and idxs[2] == localRange:upper(2) then
      SM = self._matrixFn_TR(self.dx,
                             DxxCPtr, DyyCPtr, DxyCPtr,
                             DxxLPtr, DyyLPtr, DxyLPtr,
                             DxxRPtr, DyyRPtr, DxyRPtr,
                             DxxBPtr, DyyBPtr, DxyBPtr,
                             DxxTPtr, DyyTPtr, DxyTPtr,
                             self.bcUpper[1].D,
                             self.bcUpper[1].N,
                             self.bcUpper[1].val,
                             self.bcUpper[2].D,
                             self.bcUpper[2].N,
                             self.bcUpper[2].val)
   elseif idxs[1] == localRange:lower(1) then
      SM = self._matrixFn_L(self.dx,
                            DxxCPtr, DyyCPtr, DxyCPtr,
                            DxxLPtr, DyyLPtr, DxyLPtr,
                            DxxRPtr, DyyRPtr, DxyRPtr,
                            DxxBPtr, DyyBPtr, DxyBPtr,
                            DxxTPtr, DyyTPtr, DxyTPtr,
                            self.bcLower[1].D,
                            self.bcLower[1].N,
                            self.bcLower[1].val,
                            0, 0, 0)
   elseif idxs[1] == localRange:upper(1) then
      SM = self._matrixFn_R(self.dx,
                            DxxCPtr, DyyCPtr, DxyCPtr,
                            DxxLPtr, DyyLPtr, DxyLPtr,
                            DxxRPtr, DyyRPtr, DxyRPtr,
                            DxxBPtr, DyyBPtr, DxyBPtr,
                            DxxTPtr, DyyTPtr, DxyTPtr,
                            self.bcUpper[1].D,
                            self.bcUpper[1].N,
                            self.bcUpper[1].val,
                            0, 0, 0)
   elseif idxs[2] == localRange:lower(2)then
      SM = self._matrixFn_B(self.dx,
                            DxxCPtr, DyyCPtr, DxyCPtr,
                            DxxLPtr, DyyLPtr, DxyLPtr,
                            DxxRPtr, DyyRPtr, DxyRPtr,
                            DxxBPtr, DyyBPtr, DxyBPtr,
                            DxxTPtr, DyyTPtr, DxyTPtr,
                            0, 0, 0,
                            self.bcLower[2].D,
                            self.bcLower[2].N,
                            self.bcLower[2].val)
   elseif idxs[2] == localRange:upper(2) then
      SM = self._matrixFn_T(self.dx,
                            DxxCPtr, DyyCPtr, DxyCPtr,
                            DxxLPtr, DyyLPtr, DxyLPtr,
                            DxxRPtr, DyyRPtr, DxyRPtr,
                            DxxBPtr, DyyBPtr, DxyBPtr,
                            DxxTPtr, DyyTPtr, DxyTPtr,
                            0, 0, 0,
                            self.bcUpper[2].D,
                            self.bcUpper[2].N,
                            self.bcUpper[2].val)
   else
      SM = self._matrixFn(self.dx,
                          DxxCPtr, DyyCPtr, DxyCPtr,
                          DxxLPtr, DyyLPtr, DxyLPtr,
                          DxxRPtr, DyyRPtr, DxyRPtr,
                          DxxBPtr, DyyBPtr, DxyBPtr,
                          DxxTPtr, DyyTPtr, DxyTPtr,
                          0, 0, 0, 0, 0, 0)
   end
   return SM
end

-- Build stiffness matrix
function DiscontGenPoisson:buildStiffMatrix()
   local ndim = self.ndim

   -- stencil range is a box from [-1,-1] to [1, 1] and used to index
   -- stencil elements
   local lower, upper = {}, {}
   for d = 1,ndim do
      lower[d] = -1
      upper[d] = 1
   end
   local stencilRange = Range.Range(lower, upper)
   -- we must use rowMajor indexing for stencils to be consistent with
   -- what getBlock() method returns
   local stencilIndexer = Range.makeRowMajorGenIndexer(stencilRange)

   local stiffMatrixIndexer = self.stiffMatrixIndexer
   local idxsExtRow, idxsExtCol = Lin.Vec(ndim+1), Lin.Vec(ndim+1)

   for idxs in self.Dxx:localRange():rowMajorIter() do
      local SM = self:getBlock(idxs)

      for stencilIdx in stencilRange:rowMajorIter() do
         for d = 1,ndim do
            idxsExtRow[d] = idxs[d]
            idxsExtCol[d] = idxs[d] + stencilIdx[d]
         end
         local SMij = SM[stencilIndexer(stencilIdx)]

         for k = 1,self.numBasis do
            idxsExtRow[ndim+1] = k
            local idxRow = stiffMatrixIndexer(idxsExtRow)
            for l = 1,self.basis:numBasis() do
               idxsExtCol[ndim+1] = l
               local idxCol = stiffMatrixIndexer(idxsExtCol)
               local val = SMij[k][l]
               if math.abs(val) > 1.e-14 then
                  ffiC.discontPoisson_pushTriplet(self.poisson, idxRow-1, idxCol-1, val)
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

   local stiffMatrixIndexer = self.stiffMatrixIndexer
   local srcIndexer = src:genIndexer()
   local solIndexer = sol:genIndexer()

   local idxsExt = Lin.Vec(ndim+1)
   local srcMod = Lin.Vec(self.numBasis)
   -- construct RHS vector
   for idxs in src:localRange():rowMajorIter() do
      local SM = self:getBlock(idxs)
      for k = 1,self.numBasis do
         srcMod[k] = -SM[10][k]
      end

      local srcPtr = src:get(srcIndexer(idxs))
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_pushSource(self.poisson, idx-1, srcPtr:data(), srcMod:data())
   end

   -- solve linear system
   ffiC.discontPoisson_solve(self.poisson)

   -- copy solution to field
   for idxs in src:localRange():rowMajorIter() do
      local solPtr = sol:get(solIndexer(idxs))
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_getSolution(self.poisson, idx-1, solPtr:data())
   end
end

function DiscontGenPoisson:delete()
   ffiC.delete_DiscontPoisson(self.poisson)
end

return DiscontGenPoisson
