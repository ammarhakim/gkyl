-- Gkyl ------------------------------------------------------------------------
--
-- Updater to (directly) solve Poisson equation
--
--      - Laplacian(phi) = rho
--
-- in 1D, 2D or 3D with RDG discretization.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local DataStruct = require "DataStruct"
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
local DiscontPoisson = Proto(UpdaterBase)

function DiscontPoisson:init(tbl)
   DiscontPoisson.super.init(self, tbl)

   self.grid = assert(tbl.onGrid, "Updater.DiscontPoisson: Must provide grid object using 'onGrid'")
   self.basis = assert(tbl.basis, "Updater.DiscontPoisson: Must specify basis functions to use using 'basis'")

   assert(self.grid:ndim() == self.basis:ndim(), "Dimensions of basis and grid must match")
   self.ndim = self.grid:ndim()
   self.nbasis = self.basis:numBasis()
   local polyOrder = self.basis:polyOrder()
   self.ncell = ffi.new("int[3]")
   local dx = Lin.Vec(3)    -- Limited to uniform grids for now.
   for d = 1,self.ndim do
      self.ncell[d-1] = self.grid:numCells(d)
      dx[d] = self.grid:dx(d)
   end

   local writeMatrix = xsys.pickBool(tbl.writeMatrix, false)

   -- Read the boundary conditions in
   assert(#tbl.bcLower == self.ndim, "Updater.DiscontPoisson: Must provide lower boundary conditions for all the dimesions using 'bcLower'; use empty table, '{}', for periodic BCs")
   assert(#tbl.bcUpper == self.ndim, "Updater.DiscontPoisson: Must provide upper boundary conditions for all the dimesions using 'bcUpper'; use empty table, '{}', for periodic BCs")
   local bcLower, bcUpper = {}, {}
   self.isAllPeriodic = true
   for d = 1, self.ndim do
      if self.grid:isDirPeriodic(d) then
         bcLower[d] = nil
      elseif tbl.bcLower[d].T == "D" then -- Purely Dirichlet BC
         bcLower[d] = {1, 0, tbl.bcLower[d].V}
      elseif tbl.bcLower[d].T == "N" then -- Purely Neumann BC
         bcLower[d] = {0, 1, tbl.bcLower[d].V}
      else -- Robin (mixed) BC
         bcLower[d] = {tbl.bcLower[d].D, tbl.bcLower[d].N, tbl.bcLower[d].val}
      end

      if self.grid:isDirPeriodic(d) then
         bcUpper[d] = nil
      elseif tbl.bcUpper[d].T == "D" then -- Purely Dirichlet BC
         bcUpper[d] = {1, 0, tbl.bcUpper[d].V}
      elseif tbl.bcUpper[d].T == "N" then -- Purely Neumann BC
         bcUpper[d] = {0, 1, tbl.bcUpper[d].V}
      else -- Robin (mixed) BC
         bcUpper[d] = {tbl.bcUpper[d].D, tbl.bcUpper[d].N, tbl.bcUpper[d].val}
      end

      if not self.grid:isDirPeriodic(d) then
         self.isAllPeriodic = false
      end
   end

   local basisNm = ''
   if self.ndim > 1 then
      if self.basis:id() == 'serendipity' then
         basisNm = 'Ser'
      elseif self.basis:id() == 'tensor' then
         basisNm = 'Tensor'
      end
   end

   self._first = true
   local dirs = {'x','y','z'}
   self.stencilMatrix = {}
   self.stencilMatrixLo, self.stencilMatrixUp = {}, {}
   for d = 1, self.ndim do
      local stencilMatrixFn = nil
      if polyOrder == 1 or self.ndim == 1 then 
         stencilMatrixFn = require(string.format("Updater.discontPoissonData.discontPoissonStencil%dD_%dp_%s",
                                                 self.ndim, polyOrder, dirs[d]))
      else
         stencilMatrixFn = require(string.format("Updater.discontPoissonData.discontPoisson%sStencil%dD_%dp_%s",
                                                 basisNm, self.ndim, polyOrder, dirs[d]))
      end
      self.stencilMatrix[d] = stencilMatrixFn(dx)
   end
   self.nnonzero = 0
   for k = 1,self.nbasis do
      for l = 1,self.nbasis do
         if self.stencilMatrix[1][2][k][l] ~= 0 then self.nnonzero = self.nnonzero + 1 end
         for d = 1, self.ndim do
            if self.stencilMatrix[d][1][k][l] ~= 0 then self.nnonzero = self.nnonzero + 1 end
            if self.stencilMatrix[d][3][k][l] ~= 0 then self.nnonzero = self.nnonzero + 1 end
         end
      end
   end
   for d = 1, self.ndim do
      local stencilMatrixLoFn = nil
      local stencilMatrixUpFn = nil
      if polyOrder == 1 or self.ndim == 1 then
         stencilMatrixLoFn = require(string.format("Updater.discontPoissonData.discontPoissonStencil%dD_%dp_%sLo",
                                                   self.ndim, polyOrder, dirs[d]))
         stencilMatrixUpFn = require(string.format("Updater.discontPoissonData.discontPoissonStencil%dD_%dp_%sUp",
                                                   self.ndim, polyOrder, dirs[d]))
      else
         stencilMatrixLoFn = require(string.format("Updater.discontPoissonData.discontPoisson%sStencil%dD_%dp_%sLo",
                                                   basisNm, self.ndim, polyOrder, dirs[d]))
         stencilMatrixUpFn = require(string.format("Updater.discontPoissonData.discontPoisson%sStencil%dD_%dp_%sUp",
                                                   basisNm, self.ndim, polyOrder, dirs[d]))
      end
      if self.grid:isDirPeriodic(d) then
         self.stencilMatrixLo[d] = nil
         self.stencilMatrixUp[d] = nil
      else
         self.stencilMatrixLo[d] = stencilMatrixLoFn(dx, bcLower[d][1], bcLower[d][2], bcLower[d][3])
         self.stencilMatrixUp[d] = stencilMatrixUpFn(dx, bcUpper[d][1], bcUpper[d][2], bcUpper[d][3])
      end
   end

   self.poisson = ffiC.new_DiscontPoisson(GKYL_OUT_PREFIX, self.ncell, self.ndim, self.nbasis,
                                          self.nnonzero, polyOrder, writeMatrix)

   self.dynVec = DataStruct.DynVector { numComponents = 1 }
   
   return self
end

-- Assemble the left-side matrix of the Poisson equation
function DiscontPoisson:buildStiffMatrix(phi)
   local ndim = self.ndim

   local localRange = phi:localRange()
   local lower, upper = {}, {}
   for d = 1,ndim do
      lower[d] = localRange:lower(d)
      upper[d] = localRange:upper(d)
   end
   lower[ndim+1] = 1
   upper[ndim+1] = self.nbasis
   local stiffMatrixRange = Range.Range(lower, upper)
   local stiffMatrixIndexer = Range.makeRowMajorGenIndexer(stiffMatrixRange)

   local idxsExtK, idxsExtL = Lin.Vec(ndim+1), Lin.Vec(ndim+1)
   local val, idxK, idxL
   local isLoCorner
   
   for idxs in localRange:colMajorIter() do
      isLoCorner = self.isAllPeriodic
      for d = 1,ndim do
         idxsExtK[d] = idxs[d]
         idxsExtL[d] = idxs[d]
         if idxs[d] ~= localRange:lower(d) then isLoCorner = false end
      end
      for k = 1,self.nbasis do
         idxsExtK[ndim+1] = k
         idxK = stiffMatrixIndexer(idxsExtK)
         for l = 1,self.basis:numBasis() do
            idxsExtL[ndim+1] = l
            idxL = stiffMatrixIndexer(idxsExtL)

            -- Diagonal blocks
            if isLoCorner and k == 1 and l == 1 then 
               val = 0.0
            else
               val = 0.0
               for d = 1,ndim do
                  if idxs[d] == localRange:lower(d) and not self.grid:isDirPeriodic(d) then
                     val = val + self.stencilMatrixLo[d][2][k][l]
                  elseif idxs[d] == localRange:upper(d) and not self.grid:isDirPeriodic(d) then
                     val = val + self.stencilMatrixUp[d][2][k][l]
                  else
                     val = val + self.stencilMatrix[d][2][k][l]
                  end
               end
            end
            if val ~= 0 then
               ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
            end

            -- Off-diagonal blocks
            for d = 1,ndim do
               if idxs[d] == localRange:lower(d) and not self.grid:isDirPeriodic(d) then
                  idxsExtL[d] = idxsExtL[d]+1
                  idxL = stiffMatrixIndexer(idxsExtL)
                  val = self.stencilMatrixLo[d][3][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                  end
               elseif idxs[d] == localRange:upper(d) and not self.grid:isDirPeriodic(d) then
                  idxsExtL[d] = idxs[d]-1
                  idxL = stiffMatrixIndexer(idxsExtL)
                  val = self.stencilMatrixUp[d][1][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                  end
               else -- Interior cell
                  if idxs[d] == localRange:lower(d) then
                     idxsExtL[d] = localRange:upper(d)
                  else
                     idxsExtL[d] = idxs[d]-1
                  end
                  idxL = stiffMatrixIndexer(idxsExtL)
                  val = self.stencilMatrix[d][1][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                  end

                  if idxs[d] == localRange:upper(d) then
                     idxsExtL[d] = localRange:lower(d)
                  else
                     idxsExtL[d] = idxs[d]+1
                  end
                  idxL = stiffMatrixIndexer(idxsExtL)
                  val = self.stencilMatrix[d][3][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                  end
               end
               idxsExtL[d] = idxs[d]
            end
         end
      end
   end
   ffiC.discontPoisson_constructStiffMatrix(self.poisson);
end

-- Advance method.
function DiscontPoisson:_advance(tCurr, inFld, outFld)
   local ndim = self.ndim

   local src = assert(inFld[1], "DiscontPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "DiscontPoisson.advance: Must specify an output field")

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

   local intSrcVol  = {0.0}
   -- If all directions periodic need to adjust source so that integral is 0 
   if self.isAllPeriodic then
     -- integrate source
     if self._first then
       self.calcInt = CartFieldIntegratedQuantCalc {
         onGrid = self.grid,
         basis = self.basis,
         numComponents = 1,
         quantity = "V",
       }
     end
     self.calcInt:advance(0.0, {src}, {self.dynVec})
     _, intSrcVol = self.dynVec:lastData()
   end
   local srcPeriodicMod = intSrcVol[1]/self.grid:gridVolume()*math.sqrt(2)^self.ndim

   -- Pushing source to the Eigen matrix.
   local idxsExt = Lin.Vec(ndim+1)
   local srcMod = Lin.Vec(self.nbasis)
   for idxs in localRange:colMajorIter() do
      for k = 1,self.nbasis do srcMod[k] = 0.0 end
      srcMod[1] = srcPeriodicMod
      for d = 1,ndim do
         if not self.grid:isDirPeriodic(d) then
            if idxs[d] == 1 then
               for k = 1,self.nbasis do srcMod[k] = srcMod[k] - self.stencilMatrixLo[d][1][k] end
            elseif idxs[d] == self.grid:numCells(d) then
               for k = 1,self.nbasis do srcMod[k] = srcMod[k] - self.stencilMatrixUp[d][3][k] end
            end
         end
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
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_getSolution(self.poisson, idx-1, solPtr:data())
   end

   self._first = false
end

function DiscontPoisson:delete()
   ffiC.delete_DiscontPoisson(self.poisson)
end

return DiscontPoisson
