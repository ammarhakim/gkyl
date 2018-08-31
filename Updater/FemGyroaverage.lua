-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson equation in perpendicular directions with FEM scheme
-- Perpendicular directions assumed to be first two configuration-space directions
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Proto = require "Lib.Proto"
local ffi = require "ffi"
local xsys = require "xsys"
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
  typedef struct FemGyroaverage FemGyroaverage;
  FemGyroaverage* new_FemGyroaverage(int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix);
  void delete_FemGyroaverage(FemGyroaverage* f);
  void makeGlobalStiffGy(FemGyroaverage* f, double *modifierWeight, int idx, int idy);
  void makeGyavgMatrix(FemGyroaverage* f, double *rho1, double *rho2, double *rho3, int idx);
  void finishGlobalStiffGy(FemGyroaverage* f);
  void createGlobalSrcGy(FemGyroaverage* f, double* localSrcPtr, int idx, int idy);
  void zeroGlobalSrcGy(FemGyroaverage* f);
  void allreduceGlobalSrcGy(FemGyroaverage* f, MPI_Comm comm);
  void allgatherGlobalStiffGy(FemGyroaverage* f, MPI_Comm comm);
  void getSolutionGy(FemGyroaverage* f, double* localSolPtr, int idx, int idy);
  void getNodalSolutionGy(FemGyroaverage* f, double* localSolPtr, int idx, int idy);
  void solveGy(FemGyroaverage* f);
]]
local DIRICHLET = 0
local NEUMANN = 1
local DIRICHLET_VARIABLE = 2
--local bcdata
--local mt = {}
--bcdata = ffi.metatype("bcdata_t", mt)

-- FEM Poisson solver updater object
local FemGyroaverage = Proto(UpdaterBase)

function FemGyroaverage:init(tbl)
   FemGyroaverage.super.init(self, tbl)

   -- read data from input file
   self._grid = assert(tbl.onGrid, "Updater.FemGyroaverage: Must provide grid object using 'onGrid'")
   self._phaseGrid = assert(tbl.phaseGrid, "Updater.FemGyroaverage: Must provide phaseGrid object using 'phaseGrid'")
   self._basis = assert(tbl.confBasis, "Updater.FemGyroaverage: Must specify basis functions to use using 'basis'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.FemGyroaverage: Must specify phaseBasis functions to use using 'phaseBasis'")

   assert(self._basis:id()=="serendipity", "Updater.FemGyroaverage only implemented for modal serendipity basis")

   assert(self._grid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._cdim = self._grid:ndim()
   self._ndim = self._phaseGrid:ndim()

   -- boundary conditions
   -- set C flags to indicate which directions are periodic (0-based)
   self._isDirPeriodic = ffi.new("bool[2]")
   self._isDirPeriodic[0] = self._grid:isDirPeriodic(1)
   self._isDirPeriodic[1] = self._grid:isDirPeriodic(2)

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
     if self._isDirPeriodic[dir] == false and not (self._bc[dir][0].isSet and self._bc[dir][1].isSet) then
       -- if not periodic and no BCs specified, use neumann by default 
       self._bc[dir][0] = getBcData({ T = "N", V = 0.0 })
       self._bc[dir][1] = getBcData({ T = "N", V = 0.0 })
     elseif self._isDirPeriodic[dir] then
       assert(not (self._bc[dir][0].isSet or self._bc[dir][1].isSet), "Cannot specify BCs if direction is periodic")
     end
   end

   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self._p = self._basis:polyOrder()
   self._dx = self._grid:dx(1)
   self._dy = self._grid:dx(2)

   assert(self._p == 1 or self._p == 2, "This solver only implemented for polyOrder = 1 or 2")
   assert(self._cdim == 2 or self._cdim == 3, "This solver only implemented for 2D or 3D (with no solve in 3rd dimension)")

   self._gyavg = {}
   self._first = true
   self._makeStiff = true

   -- set up constant dummy field
   self.unitWeight = DataStruct.Field {
	onGrid = self._grid,
	numComponents = self._basis:numBasis(),
	ghost = {1, 1},
   }
   local initUnit = ProjectOnBasis {
      onGrid = self._grid,
      basis = self._basis,
      evaluate = function (t,xn)
                    return 1.0
                 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,0.,{},{self.unitWeight})

   -- set up field for non-uniform and/or time-dependent modifier weight
   self.modifierWeight = DataStruct.Field {
        onGrid = self._grid,
        numComponents = self._basis:numBasis(),
        ghost = {1, 1},
   }
   -- initialize to unity
   self.modifierWeight:copy(self.unitWeight)

   -- metric-dependent coefficient fields
   self.rho1 = DataStruct.Field {
        onGrid = self._phaseGrid,
        numComponents = self._phaseBasis:numBasis(),
        ghost = {1, 1},
   }
   self.rho2 = DataStruct.Field {
        onGrid = self._phaseGrid,
        numComponents = self._phaseBasis:numBasis(),
        ghost = {1, 1},
   }
   self.rho3 = DataStruct.Field {
        onGrid = self._phaseGrid,
        numComponents = self._phaseBasis:numBasis(),
        ghost = {1, 1},
   }
   assert(tbl.rho1 and tbl.rho2 and tbl.rho3, "FemGyroaverage: must specify rho1, rho2, and rho3 coefficient fields")
   self.rho1:copy(tbl.rho1)
   self.rho2:copy(tbl.rho2)
   self.rho3:copy(tbl.rho3)

   if GKYL_HAVE_MPI then
     -- split communicators in z
     local commSet = self._grid:commSet()
     local worldComm = commSet.comm
     local nodeComm = commSet.nodeComm
     local nodeRank = Mpi.Comm_rank(nodeComm)
     local zrank = 0
     if self._cdim==3 then zrank = math.floor(nodeRank/self._grid:cuts(1)/self._grid:cuts(2)) end
     self._zcomm = Mpi.Comm_split(worldComm, zrank, nodeRank)
   end

   self.dynVec = DataStruct.DynVector { numComponents = 1 }

   return self
end

-- for testing
function FemGyroaverage:bcType(dir,side) return self._bc[dir][side].type end
function FemGyroaverage:bcValue(dir,side) return self._bc[dir][side].value end

---- advance method
function FemGyroaverage:_advance(tCurr, dt, inFld, outFld) 
   local grid = self._grid
   local basis = self._basis

   local src = assert(inFld[1], "FemGyroaverage.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemGyroaverage.advance: Must specify an output field")

   local ndim = self._ndim
   local cdim = self._cdim

   assert(#sol == self._phaseGrid:numCells(ndim), "FemGyroaverage.advance: Output field must be a length Nmu array of fields")

   -- create region that is effectively 2d and global in x-y directions
   local perpRange = self.rho1:globalRange()
   local localRange = self.rho1:localRange()
   local local_mu_lower = localRange:lower(ndim)
   local local_mu_upper = localRange:upper(ndim)
   local local_z_lower = 1
   local local_z_upper = 1
   -- remove last dimension == mu dimension from perpRange
   perpRange:shorten(ndim)
   if(cdim==3) then
     -- remove last config dimension == z dimension from perpRange
     perpRange:shorten(3)
     local_z_lower = localRange:lower(3)
     local_z_upper = localRange:upper(3)
   end

   -- create indexers and pointers for src and sol
   self.srcIndexer = src:indexer() 
   self.srcPtr = src:get(1)
   self.solIndexer = sol[1]:indexer() 
   self.solPtr = sol[1]:get(1)

   -- create indexers and pointer for other fields that belong to this object and do not change
   if self._first then 
      self.modifierWeightIndexer = self.modifierWeight:indexer() 
      self.modifierWeightPtr = self.modifierWeight:get(1)
      self.rho1Indexer = self.rho1:indexer()
      self.rho1Ptr = self.rho1:get(1)
      self.rho2Indexer = self.rho2:indexer()
      self.rho2Ptr = self.rho2:get(1)
      self.rho3Indexer = self.rho3:indexer()
      self.rho3Ptr = self.rho3:get(1)
   end

   -- loop over local z cells
   for idz=local_z_lower, local_z_upper do
     if self._first then
       self._gyavg[idz] = {}
     end
     -- loop over local mu cells
     for idmu=local_mu_lower, local_mu_upper do
       if self._first then 
         -- if first time, create gyavg C object for each z,mu cell
         self._gyavg[idz][idmu] = ffi.C.new_FemGyroaverage(self._nx, self._ny, self._cdim, self._p, 
                                             self._isDirPeriodic, self._bc, self._writeMatrix)
       end

       -- zero global source
       ffi.C.zeroGlobalSrcGy(self._gyavg[idz][idmu])

       -- create global source 
       -- globalSrc is an Eigen vector managed in C
       -- each proc allocates a full globalSrc vector
       -- loop over x and y cells locally to get local contributions to globalSrc
       for idx=localRange:lower(1),localRange:upper(1) do     
         for idy=localRange:lower(2),localRange:upper(2) do
           -- on first time, make Gyavg matrix, which varies only in x and does not vary in time
           if idy == localRange:lower(2) and self._first then
              -- note: no vpar dependence of rho's, so just use idvpar=0
              if cdim==2 then 
                self.rho1:fill(self.rho1Indexer(idx,idy,0,idmu), self.rho1Ptr) 
                self.rho2:fill(self.rho2Indexer(idx,idy,0,idmu), self.rho2Ptr) 
                self.rho3:fill(self.rho3Indexer(idx,idy,0,idmu), self.rho3Ptr) 
              else 
                self.rho1:fill(self.rho1Indexer(idx,idy,idz,0,idmu), self.rho1Ptr) 
                self.rho2:fill(self.rho2Indexer(idx,idy,idz,0,idmu), self.rho2Ptr) 
                self.rho3:fill(self.rho3Indexer(idx,idy,idz,0,idmu), self.rho3Ptr) 
              end
              ffi.C.makeGyavgMatrix(self._gyavg[idz][idmu], self.rho1Ptr:data(), self.rho2Ptr:data(), self.rho3Ptr:data(), idx-1)
           end

           if cdim==2 then 
             src:fill(self.srcIndexer(idx,idy,idmu), self.srcPtr) 
           else 
             src:fill(self.srcIndexer(idx, idy, idz, idmu), self.srcPtr) 
           end
           ffi.C.createGlobalSrcGy(self._gyavg[idz][idmu], self.srcPtr:data(), idx-1, idy-1)

           if self._makeStiff then
              if cdim==2 then 
                self.modifierWeight:fill(self.modifierWeightIndexer(idx,idy,idmu), self.modifierWeightPtr) 
              else 
                self.modifierWeight:fill(self.modifierWeightIndexer(idx, idy, idz, idmu), self.modifierWeightPtr) 
              end
              ffi.C.makeGlobalStiffGy(self._gyavg[idz][idmu], self.modifierWeightPtr:data(), idx-1, idy-1)

           end 
         end
       end

       if GKYL_HAVE_MPI and Mpi.Comm_size(self._zcomm)>1 then
         -- sum each proc's globalSrc to get final globalSrc (on each proc via allreduce)
         ffi.C.allreduceGlobalSrcGy(self._gyavg[idz][idmu], Mpi.getComm(self._zcomm))
         if self._makeStiff then
           ffi.C.allgatherGlobalStiffGy(self._gyavg[idz][idmu], Mpi.getComm(self._zcomm))
         end
       end

       if self._makeStiff then
          ffi.C.finishGlobalStiffGy(self._gyavg[idz][idmu])
       end
       -- solve 
       ffi.C.solveGy(self._gyavg[idz][idmu])

       -- remap global nodal solution to local modal solution 
       -- only need to loop over local proc region
       for idx=localRange:lower(1),localRange:upper(1) do     
         for idy=localRange:lower(2),localRange:upper(2) do
           if cdim==2 then 
             sol[idmu]:fill(self.solIndexer(idx,idy), self.solPtr) 
           else 
             sol[idmu]:fill(self.solIndexer(idx, idy, idz), self.solPtr) 
           end
           ffi.C.getSolutionGy(self._gyavg[idz][idmu], self.solPtr:data(), idx-1, idy-1)
         end
       end
     end 
   end

   self._first = false
   -- reset makeStiff flag to false, until stiffness matrix changes again
   self._makeStiff = false

   return true, GKYL_MAX_DOUBLE
end

function FemGyroaverage:delete()
  ffi.C.delete_FemGyroaverage(self._gyavg)
end

function FemGyroaverage:setLaplacianWeight(weight)
   self._hasLaplacian = true
   self.laplacianWeight:copy(weight)
   -- need to remake stiffness matrix since laplacianWeight has changed
   self._makeStiff = true
end
function FemGyroaverage:setModifierWeight(weight)
   self._hasModifier = true
   self.modifierWeight:copy(weight)
   -- need to remake stiffness matrix since modifierWeight has changed
   self._makeStiff = true
end

return FemGyroaverage
