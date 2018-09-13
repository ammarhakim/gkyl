-- Gkyl ------------------------------------------------------------------------
--
-- Updater to gyroaverage and integrate a distribution function, resulting in a gyroaveraged moment
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
local Basis = require "Basis"
local Grid = require "Grid"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local DataStruct = require "DataStruct"
local Lin = require "Lib/Linalg"
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

          int istart[8];
          int iend[8];
          int cornerstart[8];
          int cornerend[8];
      } bcdata_t;
  typedef struct FemGyroaverage FemGyroaverage;
  FemGyroaverage* new_FemGyroaverage(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix);
  void delete_FemGyroaverage(FemGyroaverage* f);
  void makeGlobalStiffGy(FemGyroaverage* f, double *modifierWeight, int idx, int idy);
  void makeGyavg0Matrix(FemGyroaverage* f, double *rho1, double *rho2, double *rho3, int idx);
  void makeGyavgMatrix(FemGyroaverage* f, double *rho1, double *rho2, double *rho3, int idx);
  void finishGlobalStiffGy(FemGyroaverage* f);
  void createGlobalSrcGy(FemGyroaverage* f, double* localSrcPtr, int idx, int idy, int ndimSrc);
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

local function isGkMomentNameGood(nm)
   if nm == "GkM0" or nm == "GkM1" or nm == "GkM2par" or nm == "GkM2perp" or nm == "GkM2" or nm == "GkM3par" or nm == "GkM3perp" or nm == "GkThreeMoments" or nm == "GkStarMoments" then
      return true
   end
   return false
end

-- FEM Poisson solver updater object
local FemGyavgGkMomentCalc = Proto(UpdaterBase)

function FemGyavgGkMomentCalc:init(tbl)
   FemGyavgGkMomentCalc.super.init(self, tbl)

   -- read data from input file
   self._grid = assert(tbl.onGrid, "Updater.FemGyavgGkMomentCalc: Must provide grid object using 'onGrid'")
   self._phaseGrid = assert(tbl.phaseGrid, "Updater.FemGyavgGkMomentCalc: Must provide phaseGrid object using 'phaseGrid'")
   self._basis = assert(tbl.confBasis, "Updater.FemGyavgGkMomentCalc: Must specify basis functions to use using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.FemGyavgGkMomentCalc: Must specify phaseBasis functions to use using 'phaseBasis'")

   assert(self._basis:id()=="serendipity", "Updater.FemGyavgGkMomentCalc only implemented for modal serendipity basis")

   assert(self._grid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._cdim = self._grid:ndim()
   self._pdim = self._phaseGrid:ndim()
   self._vdim = self._pdim - self._cdim

   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self._p = self._basis:polyOrder()
   self._dx = self._grid:dx(1)
   self._dy = self._grid:dx(2)

   assert(self._p == 1 or self._p == 2, "This solver only implemented for polyOrder = 1 or 2")
   assert(self._cdim == 2 or self._cdim == 3, "This solver only implemented for 2D or 3D (with no solve in 3rd dimension)")
   if self._p==2 then assert(self._nx>1 and self._ny>1, "Must use nx>1 and ny>1 for p=2") end

   -- create cdim+1 dimensional grid 
   local lower, upper, cells = {}, {}, {}
   for d = 1, self._cdim do
      table.insert(lower, self._grid:lower(d))
      table.insert(upper, self._grid:upper(d))
      table.insert(cells, self._grid:numCells(d))
   end
   table.insert(lower, self._phaseGrid:lower(self._pdim))
   table.insert(upper, self._phaseGrid:upper(self._pdim))
   table.insert(cells, self._phaseGrid:numCells(self._pdim))
   self._gyGrid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
   }
   -- create basis for gyGrid
   self._gyBasis = Basis.CartModalSerendipity { ndim = self._cdim+1, polyOrder = self._p }

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

   self._gyavg = {}
   self._first = true
   self._makeStiff = true

   self.f_dvpar = DataStruct.Field {
        onGrid = self._gyGrid,
        numComponents = self._gyBasis:numBasis(),
        ghost = {1, 1},
   }

   -- set up constant dummy field
   self.unitWeight = DataStruct.Field {
	onGrid = self._gyGrid,
	numComponents = self._gyBasis:numBasis(),
	ghost = {1, 1},
   }
   local initUnit = ProjectOnBasis {
      onGrid = self._gyGrid,
      basis = self._gyBasis,
      evaluate = function (t,xn)
                    return 1.0
                 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,0.,{},{self.unitWeight})

   -- set up field for non-uniform and/or time-dependent modifier weight
   self.modifierWeight = DataStruct.Field {
        onGrid = self._gyGrid,
        numComponents = self._gyBasis:numBasis(),
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
   assert(tbl.rho1 and tbl.rho2 and tbl.rho3, "FemGyavgGkMomentCalc: must specify rho1, rho2, and rho3 coefficient fields")
   self.rho1:copy(tbl.rho1)
   self.rho2:copy(tbl.rho2)
   self.rho3:copy(tbl.rho3)

   self.mass = tbl.gkfacs[1]
   self.bmag = tbl.gkfacs[2]
   self.bmagIndexer = self.bmag:indexer()
   self.bmagPtr = self.bmag:get(1)

   local id = self._phaseBasis:id()
   local mom = tbl.moment or "GkM0"
   assert(isGkMomentNameGood(mom), "FemGyavgGkMomentCalc: moment name is invalid")
   self._momCalcStep1 = MomDecl.selectGkMomCalc(mom.."_step1", id, self._cdim, self._vdim, self._p)
   self._momCalcStep2 = MomDecl.selectGkMomCalc(mom.."_step2", id, self._cdim, self._vdim, self._p)

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

   self.dxv = Lin.Vec(self._pdim) -- cell shape
   self.w = Lin.Vec(self._pdim) -- phase-space cell center

   return self
end

-- for testing
function FemGyavgGkMomentCalc:bcType(dir,side) return self._bc[dir][side].type end
function FemGyavgGkMomentCalc:bcValue(dir,side) return self._bc[dir][side].value end

---- advance method
function FemGyavgGkMomentCalc:_advance(tCurr, dt, inFld, outFld) 
   local grid = self._grid
   local basis = self._basis

   local distf = assert(inFld[1], "FemGyavgGkMomentCalc.advance: Must specify an input field")
   local gyMom = assert(outFld[1], "FemGyavgGkMomentCalc.advance: Must specify an output field")

   local pdim = self._pdim
   local cdim = self._cdim
   local ndim = self._cdim + 1 -- dimension of FEM solve. additional dimension for mu

   assert(distf:ndim()==pdim, "FemGyavgGkMomentCalc.advance: input field must have phase space dimensions.")
   assert(gyMom:ndim()==cdim, "FemGyavgGkMomentCalc.advance: output field must have dimension cdim")

   -- create region that is effectively 2d and global in x-y directions
   local perpRange = self.rho1:globalRange()
   local localRange = self.rho1:localRange()
   local local_mu_lower = localRange:lower(pdim)
   local local_mu_upper = localRange:upper(pdim)
   local local_z_lower = 1
   local local_z_upper = 1
   -- remove last dimension == mu dimension from perpRange
   perpRange:shorten(pdim)
   if(cdim==3) then
     -- remove last config dimension == z dimension from perpRange
     perpRange:shorten(3)
     local_z_lower = localRange:lower(3)
     local_z_upper = localRange:upper(3)
   end

   -- create indexers and pointers for distf and sol
   self.f_dvparIndexer = self.f_dvpar:indexer() 
   self.f_dvparPtr = self.f_dvpar:get(1)
   self.distfIndexer = distf:indexer() 
   self.distfPtr = distf:get(1)
   self.gyMomIndexer = gyMom:indexer() 
   self.gyMomPtr = gyMom:get(1)

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
         self._gyavg[idz][idmu] = ffi.C.new_FemGyroaverage(self._nx, self._ny, ndim, self._p, self._dx, self._dy,
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
           if self._first and idy==localRange:lower(2) then
              -- note: no vpar dependence of rho's, so just use idvpar=1
              if cdim==2 then 
                self.rho1:fill(self.rho1Indexer(idx,idy,1,idmu), self.rho1Ptr) 
                self.rho2:fill(self.rho2Indexer(idx,idy,1,idmu), self.rho2Ptr) 
                self.rho3:fill(self.rho3Indexer(idx,idy,1,idmu), self.rho3Ptr) 
              else 
                self.rho1:fill(self.rho1Indexer(idx,idy,idz,1,idmu), self.rho1Ptr) 
                self.rho2:fill(self.rho2Indexer(idx,idy,idz,1,idmu), self.rho2Ptr) 
                self.rho3:fill(self.rho3Indexer(idx,idy,idz,1,idmu), self.rho3Ptr) 
              end
              ffi.C.makeGyavgMatrix(self._gyavg[idz][idmu], self.rho1Ptr:data(), self.rho2Ptr:data(), self.rho3Ptr:data(), idx-1) 
           end

           -- set pointer for src = int dvpar distf
           if cdim==2 then 
             self.f_dvpar:fill(self.f_dvparIndexer(idx,idy,idmu), self.f_dvparPtr)
             self.bmag:fill(self.bmagIndexer(idx,idy), self.bmagPtr) 
           else 
             self.f_dvpar:fill(self.f_dvparIndexer(idx,idy,idz,idmu), self.f_dvparPtr)
             self.bmag:fill(self.bmagIndexer(idx,idy,idz), self.bmagPtr) 
           end
           -- loop over vpar to integrate in vpar
           for idv=localRange:lower(cdim+1),localRange:upper(cdim+1) do
             if cdim==2 then 
               distf:fill(self.distfIndexer(idx,idy,idv,idmu), self.distfPtr) 
               self._phaseGrid:setIndexes(idx,idy,idv,idmu)
             else 
               distf:fill(self.distfIndexer(idx,idy,idz,idv,idmu), self.distfPtr) 
               self._phaseGrid:setIndexes(idx,idy,idz,idv,idmu)
             end
             self._phaseGrid:cellCenter(self.w)
             for d = 1, self._pdim do self.dxv[d] = self._phaseGrid:dx(d) end
             self._momCalcStep1(self.w:data(), self.dxv:data(), self.mass, self.bmagPtr:data(), self.distfPtr:data(), self.f_dvparPtr:data())
           end

           ffi.C.createGlobalSrcGy(self._gyavg[idz][idmu], self.f_dvparPtr:data(), idx-1, idy-1, self.f_dvpar:ndim())

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
             self.f_dvpar:fill(self.f_dvparIndexer(idx,idy,idmu), self.f_dvparPtr) 
             gyMom:fill(self.gyMomIndexer(idx, idy), self.gyMomPtr) 
             self.bmag:fill(self.bmagIndexer(idx,idy), self.bmagPtr) 
             self._phaseGrid:setIndexes(idx,idy,1,idmu)
           else
             self.f_dvpar:fill(self.f_dvparIndexer(idx, idy, idz, idmu), self.f_dvparPtr) 
             gyMom:fill(self.gyMomIndexer(idx, idy, idz), self.gyMomPtr) 
             self.bmag:fill(self.bmagIndexer(idx,idy,idz), self.bmagPtr) 
             self._phaseGrid:setIndexes(idx,idy,idz,1,idmu)
           end
           ffi.C.getSolutionGy(self._gyavg[idz][idmu], self.f_dvparPtr:data(), idx-1, idy-1)

           -- integrate in mu
           self._phaseGrid:cellCenter(self.w)
           for d = 1, self._pdim do self.dxv[d] = self._phaseGrid:dx(d) end
           self._momCalcStep2(self.w:data(), self.dxv:data(), self.mass, self.bmagPtr:data(), self.f_dvparPtr:data(), self.gyMomPtr:data())
         end
       end
     end 
   end

   self._first = false
   -- reset makeStiff flag to false, until stiffness matrix changes again
   self._makeStiff = false

   return true, GKYL_MAX_DOUBLE
end

function FemGyavgGkMomentCalc:delete()
  ffi.C.delete_FemGyavgGkMomentCalc(self._gyavg)
end

function FemGyavgGkMomentCalc:setLaplacianWeight(weight)
   self._hasLaplacian = true
   self.laplacianWeight:copy(weight)
   -- need to remake stiffness matrix since laplacianWeight has changed
   self._makeStiff = true
end
function FemGyavgGkMomentCalc:setModifierWeight(weight)
   self._hasModifier = true
   self.modifierWeight:copy(weight)
   -- need to remake stiffness matrix since modifierWeight has changed
   self._makeStiff = true
end

return FemGyavgGkMomentCalc
