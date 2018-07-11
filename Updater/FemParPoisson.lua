-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson equation in parallel direction with FEM scheme
-- Parallel direction assumed to be last configuration-space direction (z)
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
local xsys = require "xsys"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local DataStruct = require "DataStruct"
local Mpi
if GKYL_HAVE_MPI then Mpi = require "Comm.Mpi" end

ffi.cdef[[
/** Structure to store BC data. */
      typedef struct 
      {
/** Flag to indicate if Bc was set */
          bool isSet;
/** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
          unsigned type;
/** Value to apply */
          double value;
          int istart; 
          int iend;
      } bcdataPar_t;
  typedef struct FemParPoisson FemParPoisson;
  FemParPoisson* new_FemParPoisson(int nz, int ndim, int polyOrder, double dz, bool periodicFlg, bcdataPar_t bc[2], bool writeMatrix, double laplacianWeight, double modifierConstant);
  void delete_FemParPoisson(FemParPoisson* f);
  void createParGlobalSrc(FemParPoisson* f, double* localSrcPtr, int idz, double intSrcVol);
  void allreduceParGlobalSrc(FemParPoisson* f, MPI_Comm comm);
  void zeroParGlobalSrc(FemParPoisson* f);
  void solvePar(FemParPoisson* f);
  void getSolutionPar(FemParPoisson* f, double* ptr, int idz);
  void getNodalSolutionPar(FemParPoisson* f, double* ptr, int idz);
]]
local DIRICHLET = 0
local NEUMANN = 1
local DIRICHLET_VARIABLE = 2

-- FEM Poisson solver updater object
local FemParPoisson = Proto(UpdaterBase)

function FemParPoisson:init(tbl)
   FemParPoisson.super.init(self, tbl)

   -- read data from input file
   self._grid = assert(tbl.onGrid, "Updater.FemParPoisson: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemParPoisson: Must specify basis functions to use using 'basis'")

   assert(self._basis:id()=="serendipity", "Updater.FemParPoisson only implemented for modal serendipity basis")

   assert(self._grid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._grid:ndim()
   -- solve direction is parallel (z) direction, which is always assumed to be last config space direction
   self._zdir = self._ndim

   -- boundary conditions
   -- extract periodic directions
   self._periodic = false
   if tbl.periodicDirs then 
      for i, d in ipairs(tbl.periodicDirs) do
         if d == self._zdir then self._periodic = true end
      end
   end

   self._writeMatrix = xsys.pickBool(tbl.writeStiffnessMatrix, false)
  
   local function getBcData(tbl)
     local bc = ffi.new("bcdataPar_t")
     if tbl.T == "D" then bc.type = DIRICHLET
     elseif tbl.T == "N" then bc.type = NEUMANN
     elseif tbl.T == "D_VAR" then bc.type = DIRICHLET_VARIABLE
     else assert(false, "Boundary condition type must be specified by one of 'D', 'N', or 'D_VAR'")
     end
     bc.value = tbl.V
     bc.isSet = true
     return bc
   end

   self._bc = ffi.new("bcdataPar_t[2]")
 
   if tbl.bcBack then
     self._bc[0] = getBcData(tbl.bcBack)
   end
   if tbl.bcFront then
     self._bc[1] = getBcData(tbl.bcFront)
   end

   -- make sure BCs are specified consistently
   if self._periodic == false and not (self._bc[0].isSet and self._bc[1].isSet) then
     -- if not periodic, use neumann by default (usually doing discont-to-cont projection)
     self._bc[0] = getBcData({ T = "N", V = 0.0 })
     self._bc[1] = getBcData({ T = "N", V = 0.0 })
   elseif self._periodic then
     assert(not (self._bc[0].isSet or self._bc[1].isSet), "Cannot specify BCs if direction is periodic")
   end

   self._modifierConstant=0.0
   if tbl.modifierConstant then
     self._modifierConstant = tbl.modifierConstant
   end

   self._laplacianWeight=1.0
   if tbl.laplacianWeight then
     self._laplacianWeight = tbl.laplacianWeight
   end

   self._adjustSource = false
   if self._periodic and self._modifierConstant == 0.0 then self._adjustSource = true end

   self._nz = self._grid:numCells(self._zdir)
   self._p = self._basis:polyOrder()
   self._dz = self._grid:dx(self._zdir)

   assert(self._p == 1 or self._p == 2, "This solver only implemented for polyOrder = 1 or 2")
   assert(self._ndim <= 3, "This solver only implemented for 1D, 2D or 3D (with a solve only in last dimension)")

   self._poisson = ffi.C.new_FemParPoisson(self._nz, self._ndim, self._p, 
                                            self._dz, self._periodic, 
                                            self._bc, self._writeMatrix,
                                            self._laplacianWeight, self._modifierConstant)

   if GKYL_HAVE_MPI then
     -- split communicators in x-y
     local commSet = self._grid:commSet()
     local worldComm = commSet.comm
     local nodeComm = commSet.nodeComm
     local nodeRank = Mpi.Comm_rank(nodeComm)
     local xyrank = 0
     if self._ndim>1 then xyrank = nodeRank%(self._grid:cuts(1)*self._grid:cuts(2)) end
     self._xycomm = Mpi.Comm_split(worldComm, xyrank, nodeRank)
   end

   self.dynVec = DataStruct.DynVector { numComponents = 1 }

   return self
end

-- for testing
function FemParPoisson:bcType(side) return self._bc[side].type end
function FemParPoisson:bcValue(side) return self._bc[side].value end

---- advance method
function FemParPoisson:_advance(tCurr, dt, inFld, outFld) 
   local grid = self._grid
   local basis = self._basis

   local src = assert(inFld[1], "FemParPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemParPoisson.advance: Must specify an output field")

   local ndim = self._ndim

   -- create region that is effectively 1d and global in z directions
   local parRange = src:globalRange()
   local localRange = src:localRange()
   local local_xy_lower = {1, 1}
   local local_xy_upper = {1, 1}
   for d = 1, self._ndim-1 do
      parRange:shorten(d)
      local_xy_lower[d] = localRange:lower(d)
      local_xy_upper[d] = localRange:upper(d)
   end

   -- create indexers and pointers for src and sol
   local srcIndexer = src:indexer() 
   local srcPtr = src:get(1)
   local solIndexer = sol:indexer() 
   local solPtr = sol:get(1)

   local intSrcVol = {0.0}
   -- if all directions periodic need to adjust source so that integral is 0 
   if self._adjustSource then
     -- integrate source
     local calcInt = CartFieldIntegratedQuantCalc {
       onGrid = grid,
       basis = basis,
       numComponents = 1,
       quantity = "V",
     }
     calcInt:advance(0.0, 0.0, {src}, {self.dynVec})
     _, intSrcVol = self.dynVec:lastData()
   end

   -- loop over local x-y cells
   for idx=local_xy_lower[1], local_xy_upper[1] do
     for idy=local_xy_lower[2], local_xy_upper[2] do
       -- zero global source
       ffi.C.zeroParGlobalSrc(self._poisson)

       -- create global source 
       -- globalSrc is an Eigen vector managed in C
       -- each proc allocates a full globalSrc vector
       -- loop over z cells locally to get local contributions to globalSrc
       for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
         if ndim==1 then 
           src:fill(srcIndexer(idz), srcPtr)
         elseif ndim==2 then
           src:fill(srcIndexer(idx, idz), srcPtr)
         else
           src:fill(srcIndexer(idx, idy, idz), srcPtr)
         end
         ffi.C.createParGlobalSrc(self._poisson, srcPtr:data(), idz-1, intSrcVol[1]/grid:gridVolume()*math.sqrt(2)^ndim)
       end

       if GKYL_HAVE_MPI and Mpi.Comm_size(self._xycomm)>1 then
         -- sum each proc's globalSrc to get final globalSrc (on each proc via allreduce)
         ffi.C.allreduceParGlobalSrc(self._poisson, Mpi.getComm(self._xycomm))
       end

       -- solve 
       ffi.C.solvePar(self._poisson)
  
       -- remap global nodal solution to local modal solution 
       -- only need to loop over local proc region
       for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
         if ndim==1 then 
           sol:fill(solIndexer(idz), solPtr)
         elseif ndim==2 then
           sol:fill(solIndexer(idx, idz), solPtr)
         else
           sol:fill(solIndexer(idx, idy, idz), solPtr)
         end
         ffi.C.getSolutionPar(self._poisson, solPtr:data(), idz-1)
       end
     end
   end

   return true, GKYL_MAX_DOUBLE
end

function FemParPoisson:delete()
  ffi.C.delete_FemParPoisson(self._poisson)
end

return FemParPoisson
