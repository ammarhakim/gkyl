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

          int istart[3];
          int iend[3];
          int cornerstart[3];
          int cornerend[3];
      } bcdata_t;
  typedef struct FemPerpPoisson FemPerpPoisson;
  FemPerpPoisson* new_FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix, double laplacianWeight, double modifierConstant);
  void delete_FemPerpPoisson(FemPerpPoisson* f);
  void createGlobalSrc(FemPerpPoisson* f, double* localSrcPtr, int idx, int idy, double intSrcVol);
  void allreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm);
  void zeroGlobalSrc(FemPerpPoisson* f);
  void solve(FemPerpPoisson* f);
  void getSolution(FemPerpPoisson* f, double* ptr, int idx, int idy);
  void getNodalSolution(FemPerpPoisson* f, double* ptr, int idx, int idy);
]]
local DIRICHLET = 0
local NEUMANN = 1
local DIRICHLET_VARIABLE = 2
--local bcdata
--local mt = {}
--bcdata = ffi.metatype("bcdata_t", mt)

-- FEM Poisson solver updater object
local FemPerpPoisson = Proto(UpdaterBase)

function FemPerpPoisson:init(tbl)
   FemPerpPoisson.super.init(self, tbl)

   -- read data from input file
   self._grid = assert(tbl.onGrid, "Updater.FemPerpPoisson: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemPerpPoisson: Must specify basis functions to use using 'basis'")

   assert(self._basis:id()=="serendipity", "Updater.FemPerpPoisson only implemented for modal serendipity basis")

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

   self._writeMatrix = xsys.pickBool(tbl.writeStiffnessMatrix, false)
  
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
       assert(self._bc[dir][0].isSet and self._bc[dir][1].isSet, "Must specify non-periodic BCs on each side")
     else
       assert(not (self._bc[dir][0].isSet or self._bc[dir][1].isSet), "Cannot specify BCs if direction is periodic")
     end
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
   if self._allPeriodic and self._modifierConstant == 0.0 then self._adjustSource = true end

   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self._p = self._basis:polyOrder()
   self._dx = self._grid:dx(1)
   self._dy = self._grid:dx(2)

   assert(self._p == 1 or self._p == 2, "This solver only implemented for polyOrder = 1 or 2")
   assert(self._ndim == 2 or self._ndim == 3, "This solver only implemented for 2D or 3D (with no solve in 3rd dimension)")

   self._poisson = ffi.C.new_FemPerpPoisson(self._nx, self._ny, self._ndim, self._p, 
                                            self._dx, self._dy, self._isDirPeriodic, 
                                            self._bc, self._writeMatrix,
                                            self._laplacianWeight, self._modifierConstant)

   if GKYL_HAVE_MPI then
     -- split communicators in z
     local commSet = self._grid:commSet()
     local worldComm = commSet.comm
     local nodeComm = commSet.nodeComm
     local nodeRank = Mpi.Comm_rank(nodeComm)
     local zrank = 0
     if self._ndim==3 then zrank = math.floor(nodeRank/self._grid:cuts(1)/self._grid:cuts(2)) end
     self._zcomm = Mpi.Comm_split(worldComm, zrank, nodeRank)
   end

   self.zContinuous = xsys.pickBool(tbl.zContinuous, false)
   self.zDiscontToCont = nil
   if self._ndim == 3 and self.zContinuous then
     self.zDisContToCont = Updater.FemParPoisson {
       onGrid = self._grid,
       basis = self._basis,
       bcBack = { T = "N", V = 0.0 },
       bcFront = { T = "N", V = 0.0 },
       laplacianWeight = 0.0,
       modifierConstant = 1.0,
     }
   end

   return self
end

-- for testing
function FemPerpPoisson:bcType(dir,side) return self._bc[dir][side].type end
function FemPerpPoisson:bcValue(dir,side) return self._bc[dir][side].value end

---- advance method
function FemPerpPoisson:_advance(tCurr, dt, inFld, outFld) 
   local grid = self._grid
   local basis = self._basis

   local src = assert(inFld[1], "FemPerpPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemPerpPoisson.advance: Must specify an output field")

   local ndim = self._ndim

   -- create region that is effectively 2d and global in x-y directions
   local perpRange = src:globalRange()
   local localRange = src:localRange()
   local local_z_lower = 1
   local local_z_upper = 1
   if(ndim==3) then
     perpRange:shorten(3)
     local_z_lower = localRange:lower(3)
     local_z_upper = localRange:upper(3)
   end

   -- create indexers and pointers for src and sol
   local srcIndexer = src:indexer() 
   local srcPtr = src:get(1)
   local solIndexer = sol:indexer() 
   local solPtr = sol:get(1)

   local intSrcVol = {0.0}
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   -- if all directions periodic need to adjust source so that integral is 0 
   if self._adjustSource then
     -- integrate source
     local calcInt = CartFieldIntegratedQuantCalc {
       onGrid = grid,
       basis = basis,
       numComponents = 1,
     }
     calcInt:advance(0.0, 0.0, {src}, {dynVec})
     _, intSrcVol = dynVec:lastData()
   end

   -- loop over local z cells
   for idz=local_z_lower, local_z_upper do
     -- zero global source
     ffi.C.zeroGlobalSrc(self._poisson)

     -- create global source 
     -- globalSrc is an Eigen vector managed in C
     -- each proc allocates a full globalSrc vector
     -- loop over x and y cells locally to get local contributions to globalSrc
     for idx=localRange:lower(1),localRange:upper(1) do     
       for idy=localRange:lower(2),localRange:upper(2) do
         if ndim==2 then 
           src:fill(srcIndexer(idx,idy), srcPtr) 
         else 
           src:fill(srcIndexer(idx, idy, idz), srcPtr) 
         end
         ffi.C.createGlobalSrc(self._poisson, srcPtr:data(), idx-1, idy-1, intSrcVol[1]/grid:gridVolume())
       end
     end

     if GKYL_HAVE_MPI and Mpi.Comm_size(self._zcomm)>1 then
       -- sum each proc's globalSrc to get final globalSrc (on each proc via allreduce)
       ffi.C.allreduceGlobalSrc(self._poisson, Mpi.getComm(self._zcomm))
     end

     -- solve 
     ffi.C.solve(self._poisson)

     -- remap global nodal solution to local modal solution 
     -- only need to loop over local proc region
     for idx=localRange:lower(1),localRange:upper(1) do     
       for idy=localRange:lower(2),localRange:upper(2) do
         if ndim==2 then 
           sol:fill(solIndexer(idx,idy), solPtr) 
         else 
           sol:fill(solIndexer(idx, idy, idz), solPtr) 
         end
         ffi.C.getSolution(self._poisson, solPtr:data(), idx-1, idy-1)
       end
     end
   end 

   -- optionally make continuous in z
   if self.zDiscontToCont then self.zDiscontToCont:advance(tCurr, dt, {sol}, {sol}) end

   return true, GKYL_MAX_DOUBLE
end

function FemPerpPoisson:delete()
  ffi.C.delete_FemPerpPoisson(self._poisson)
end

return FemPerpPoisson
