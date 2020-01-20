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
local Mpi
-- if GKYL_HAVE_MPI then Mpi = require "Comm.Mpi" end

ffi.cdef[[
  /** Structure to store BC data. */
  typedef struct {
    /** Flag to indicate if Bc was set */
    bool isSet;
    /** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
    unsigned type;
    /** Value to apply */
    double value;

    int istart[8];
    int iend[8];
    int cornerstart[8];
    int cornerend[8];
  } bcdata_t;
  typedef struct DiscontPoisson DiscontPoisson;
  DiscontPoisson* new_DiscontPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2]);
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
local bcdata
local mt = {}
bcdata = ffi.metatype("bcdata_t", mt)

-- FEM Poisson solver updater object
local DiscontPoisson = Proto(UpdaterBase)

function DiscontPoisson:init(tbl)
   DiscontPoisson.super.init(self, tbl)

   -- read data from input file
   self._grid = assert(tbl.onGrid, "Updater.DiscontPoisson: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.DiscontPoisson: Must specify basis functions to use using 'basis'")

   assert(self._basis:id() == "serendipity", "Updater.DiscontPoisson only implemented for modal serendipity basis")

   assert(self._grid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._grid:ndim()

   -- boundary conditions
   -- extract periodic directions
   local periodicDirs = {}
   if tbl.periodicDirs then
      for i, d in ipairs(tbl.periodicDirs) do
         if d<1 or d>self._ndim then
            assert(false, "Directions in periodicDirs table should be 1 (for X), or 2 (for Y)")
         end
         periodicDirs[i] = d
      end
   end
   -- set C flags to indicate which directions are periodic (0-based)
   self._isDirPeriodic = ffi.new("bool[2]")
   self._isDirPeriodic[0] = false
   self._isDirPeriodic[1] = false
   for _, d in ipairs(periodicDirs) do self._isDirPeriodic[d-1] = true end

   -- set flag to indicate all directions are periodic 
   self._allPeriodic = true
   for d = 0,1 do
     if not self._isDirPeriodic[d] then 
       self._allPeriodic = false 
     end
   end

   --self._writeMatrix = xsys.pickBool(tbl.writeStiffnessMatrix, false)
  
   local function getBcData(tbl)
     local bc = ffi.new("bcdata_t")
     if tbl.T == "D" then bc.type = DIRICHLET
     elseif tbl.T == "N" then bc.type = NEUMANN
     elseif tbl.T == "D_VAR" then bc.type = DIRICHLET_VARIABLE
     else assert(false, "Boundary condition type must be specified by one of 'D', 'N', or 'D_VAR'")
     end
     bc.value = tbl.V
     bc.isSet = true
     return bc
   end

   self._bc = ffi.new("bcdata_t[2][2]")
 
   if tbl.bcLeft then
     self._bc[0][0] = getBcData(tbl.bcLeft)
   end
   if tbl.bcRight then
     self._bc[0][1] = getBcData(tbl.bcRight)
   end
   if tbl.bcBottom then
     self._bc[1][0] = getBcData(tbl.bcBottom)
   end
   if tbl.bcTop then
     self._bc[1][1] = getBcData(tbl.bcTop)
   end

   -- make sure BCs are specified consistently
   for dir=0,1 do
     if self._isDirPeriodic[dir] == false then
       assert(self._bc[dir][0].isSet and self._bc[dir][1].isSet, "Must specify non-periodic BCs on each side (dir " .. dir .. ")")
     else
       assert(not (self._bc[dir][0].isSet or self._bc[dir][1].isSet), "Cannot specify BCs if direction is periodic")
     end
   end

   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self._p = self._basis:polyOrder()
   self._dx = self._grid:dx(1)
   self._dy = self._grid:dx(2)

   assert(self._p == 1 or self._p == 2, "This solver only implemented for polyOrder = 1 or 2")
   assert(self._ndim == 2 or self._ndim == 3, "This solver only implemented for 2D or 3D")
   if self._p == 2 then assert(self._nx>1 and self._ny>1, "Must use nx>1 and ny>1 for p=2") end

   self._first = true
   self._makeStiff = true

   self._poisson = ffiC.new_DiscontPoisson(self._nx, self._ny, self._ndim, self._p, 
                                           self._dx, self._dy, self._isDirPeriodic,
                                           self._bc)

   return self
end

-- for testing
function DiscontPoisson:bcType(dir,side) return self._bc[dir][side].type end
function DiscontPoisson:bcValue(dir,side) return self._bc[dir][side].value end

---- advance method
function DiscontPoisson:_advance(tCurr, inFld, outFld) 
   local grid = self._grid
   local basis = self._basis

   local src = assert(inFld[1], "DiscontPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "DiscontPoisson.advance: Must specify an output field")

   local ndim = self._ndim

   local stencilMatrix = require(string.format("Updater.discontPoissonData.discontPoissonStencil%dx%dp", self._ndim, self._p))
   -- create region that is effectively 2d and global in x-y directions
   local globalRange = src:globalRange()
   local localRange = src:localRange()
   local stiffMatrixRange = Range.Range({1,1,1}, {localRange:upper(1), localRange:upper(2), basis:numBasis()})
   local stiffMatrixIndexer = Range.makeRowMajorIndexer(stiffMatrixRange)

   -- create indexers and pointers for src and sol
   if self._first then 
      self.srcIndexer = src:genIndexer() 
      self.solIndexer = sol:genIndexer() 
   end

   if self._first then
      for idxs in localRange:colMajorIter() do
         for k = 1,basis:numBasis() do
            for l = 1,basis:numBasis() do
               local idxK = stiffMatrixIndexer(idxs[1], idxs[2], k)
               local idxL = stiffMatrixIndexer(idxs[1], idxs[2], l)
               local val = stencilMatrix[1][k][l]
               if val ~= 0 then
                  ffiC.discontPoisson_pushTriplet(self._poisson, idxK-1, idxL-1, val)
               end

               if idxs[1] > localRange:lower(1) then
                  idxL = stiffMatrixIndexer(idxs[1]-1, idxs[2], l)
                  val = stencilMatrix[2][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self._poisson, idxK-1, idxL-1, val)
                  end
               end

               if idxs[1] < localRange:upper(1) then
                  idxL = stiffMatrixIndexer(idxs[1]+1, idxs[2], l)
                  val = stencilMatrix[3][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self._poisson, idxK-1, idxL-1, val)
                  end
               end

               if idxs[2] > localRange:lower(2) then
                  idxL = stiffMatrixIndexer(idxs[1], idxs[2]-1, l)
                  val = stencilMatrix[4][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self._poisson, idxK-1, idxL-1, val)
                  end
               end
               
               if idxs[2] < localRange:upper(2) then
                  idxL = stiffMatrixIndexer(idxs[1], idxs[2]+1, l)
                  val = stencilMatrix[5][k][l]
                  if val ~= 0 then
                     ffiC.discontPoisson_pushTriplet(self._poisson, idxK-1, idxL-1, val)
                  end
               end
            end
         end
      end
      ffiC.discontPoisson_constructStiffMatrix(self._poisson);
   end

   for idxs in localRange:colMajorIter() do
      local srcPtr = src:get(self.srcIndexer(idxs))
      local idx = stiffMatrixIndexer(idxs[1], idxs[2], 1)
      ffiC.discontPoisson_pushSource(self._poisson, idx-1, srcPtr:data())
   end

   ffiC.discontPoisson_solve(self._poisson)

   for idxs in localRange:colMajorIter() do
      local solPtr = sol:get(self.solIndexer(idxs))
      local idx = stiffMatrixIndexer(idxs[1], idxs[2], 1)
      ffiC.discontPoisson_getSolution(self._poisson, idx-1, solPtr:data())
   end

   self._first = false
end

function DiscontPoisson:delete()
  ffiC.delete_DiscontPoisson(self._poisson)
end

return DiscontPoisson
