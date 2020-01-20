-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson equation in perpendicular directions with FEM scheme
-- Perpendicular directions assumed to be first two configuration-space directions
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local DataStruct = require "DataStruct"

ffi.cdef[[
  typedef struct DiscontPoisson DiscontPoisson;
  DiscontPoisson* new_DiscontPoisson(int ncells[3], int ndim, int nbasis, int nnonzero, int polyOrder, double dx[3]);
  void delete_DiscontPoisson(DiscontPoisson* f);

  void discontPoisson_pushTriplet(DiscontPoisson* f, int idxK, int idxL, double val); 
  void discontPoisson_constructStiffMatrix(DiscontPoisson* f);
  void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src);
  void discontPoisson_solve(DiscontPoisson* f);
  void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol);
]]
local DIRICHLET = 0
local NEUMANN = 1
local DIRICHLET_VARIABLE = 2

-- FEM Poisson solver updater object
local DiscontPoisson = Proto(UpdaterBase)

function DiscontPoisson:init(tbl)
   DiscontPoisson.super.init(self, tbl)

   -- read data from input file
   self.grid = assert(tbl.onGrid, "Updater.DiscontPoisson: Must provide grid object using 'onGrid'")
   self.basis = assert(tbl.basis, "Updater.DiscontPoisson: Must specify basis functions to use using 'basis'")

   assert(self.basis:id() == "serendipity", "Updater.DiscontPoisson only implemented for modal serendipity basis")

   assert(self.grid:ndim() == self.basis:ndim(), "Dimensions of basis and grid must match")
   self.ndim = self.grid:ndim()
   self.nbasis = self.basis:numBasis()
   self.polyOrder = self.basis:polyOrder()
   self.ncell = ffi.new("int[3]")
   self.dx = ffi.new("double[3]")
   for d = 1,self.ndim do
      self.ncell[d-1] = self.grid:numCells(d)
      self.dx[d-1] = self.grid:dx(d)
   end

   self._first = true

   self.stencilMatrix = require(string.format("Updater.discontPoissonData.discontPoissonStencil%dx%dp", self.ndim, self.polyOrder))
   self.nnonzero = 0
   for d = 1,(1+2*self.ndim) do
      for k = 1,self.nbasis do
         for l = 1,self.nbasis do
            if self.stencilMatrix[d][k][l] ~= 0 then self.nnonzero = self.nnonzero + 1 end
         end
      end
   end
   print(self.nnonzero)
      
   self.poisson = ffiC.new_DiscontPoisson(self.ncell, self.ndim, self.nbasis,
                                          self.nnonzero, self.polyOrder, self.dx)

   return self
end

---- advance method
function DiscontPoisson:_advance(tCurr, inFld, outFld) 
   local grid = self.grid
   local basis = self.basis
   local ndim = self.ndim

   local src = assert(inFld[1], "DiscontPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "DiscontPoisson.advance: Must specify an output field")

   local globalRange = src:globalRange()
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

   -- construct the stiffness matrix using Eigen
   if self._first then 
      self.srcIndexer = src:genIndexer() 
      self.solIndexer = sol:genIndexer() 

      for idxs in localRange:colMajorIter() do
         local idxsExtK, idxsExtL = Lin.Vec(ndim+1), Lin.Vec(ndim+1)
         for d = 1,ndim do
            idxsExtK[d] = idxs[d]
            idxsExtL[d] = idxs[d]
         end
         for k = 1,self.nbasis do
            idxsExtK[ndim+1] = k
            local idxK = stiffMatrixIndexer(idxsExtK)               
            for l = 1,basis:numBasis() do
               idxsExtL[ndim+1] = l
               local idxL = stiffMatrixIndexer(idxsExtL)
               
               local val = self.stencilMatrix[1][k][l]
               if val ~= 0 then
                  ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
               end

               local cnt = 2
               for d = 1,ndim do
                  if idxs[d] > localRange:lower(d) then
                     idxsExtL[d] = idxsExtL[d]-1
                     idxL = stiffMatrixIndexer(idxsExtL)
                     val = self.stencilMatrix[cnt][k][l]
                     if val ~= 0 then
                        ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                     end
                     idxsExtL[d] = idxsExtL[d]+1
                  end
                  cnt = cnt + 1

                  if idxs[d] < localRange:upper(d) then
                     idxsExtL[d] = idxsExtL[d]+1
                     idxL = stiffMatrixIndexer(idxsExtL)
                     val = self.stencilMatrix[cnt][k][l]
                     if val ~= 0 then
                        ffiC.discontPoisson_pushTriplet(self.poisson, idxK-1, idxL-1, val)
                     end
                     idxsExtL[d] = idxsExtL[d]-1
                  end
                  cnt = cnt + 1
               end
            end
         end
      end
      ffiC.discontPoisson_constructStiffMatrix(self.poisson);
   end

   for idxs in localRange:colMajorIter() do
      local srcPtr = src:get(self.srcIndexer(idxs))
      local idxsExt = Lin.Vec(ndim+1)
      for d = 1,ndim do idxsExt[d] = idxs[d] end
      idxsExt[ndim+1] = 1
      local idx = stiffMatrixIndexer(idxsExt)
      ffiC.discontPoisson_pushSource(self.poisson, idx-1, srcPtr:data())
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

function DiscontPoisson:delete()
  ffiC.delete_DiscontPoisson(self.poisson)
end

return DiscontPoisson
