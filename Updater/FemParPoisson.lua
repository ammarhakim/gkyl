-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson equation in parallel direction with FEM scheme
-- Parallel direction assumed to be last configuration-space direction (z)
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc          = require "Lib.Alloc"
local UpdaterBase    = require "Updater.Base"
local Lin            = require "Lib.Linalg"
local lume           = require "Lib.lume"
local Proto          = require "Lib.Proto"
local Range          = require "Lib.Range"
local ffi            = require "ffi"
local ffiC           = ffi.C
local xsys           = require "xsys"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local DataStruct     = require "DataStruct"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
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
  FemParPoisson* new_FemParPoisson(int nz, int ndim, int polyOrder, double dz, bool periodicFlg, bcdataPar_t bc[2], bool writeMatrix);
  void delete_FemParPoisson(FemParPoisson* f);
  void makeParGlobalStiff(FemParPoisson* f, double *laplacianWeight, double *modifierWeight, int idz);
  void finishParGlobalStiff(FemParPoisson* f);
  void createParGlobalSrc(FemParPoisson* f, double* localSrcPtr, int idz, double intSrcVol);
  void allreduceParGlobalSrc(FemParPoisson* f, MPI_Comm comm);
  void allgatherParGlobalStiff(FemParPoisson* f, MPI_Comm comm);
  void zeroParGlobalSrc(FemParPoisson* f);
  void solvePar(FemParPoisson* f);
  void getSolutionPar(FemParPoisson* f, double* ptr, int idz);
  void getNodalSolutionPar(FemParPoisson* f, double* ptr, int idz);
]]
local DIRICHLET = 0
local NEUMANN   = 1
local DIRICHLET_VARIABLE = 2

-- FEM Poisson solver updater object.
local FemParPoisson = Proto(UpdaterBase)

function FemParPoisson:init(tbl)
   FemParPoisson.super.init(self, tbl)

   -- Read data from input file.
   self._grid  = assert(tbl.onGrid, "Updater.FemParPoisson: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemParPoisson: Must specify basis functions to use using 'basis'")

   self._ndim = self._grid:ndim()

   assert(self._basis:id()=="serendipity", "Updater.FemParPoisson: only implemented for modal serendipity basis")
   assert(self._basis:polyOrder()==1 or self._basis:polyOrder()==2, "Updater.FemParPoisson: only implemented for polyOrder = 1 or 2")
   assert(self._ndim == self._basis:ndim(), "Updater.FemParPoisson: Dimensions of basis and grid must match")
   assert(self._ndim==1 or self._ndim==3, "Updater.FemParPoisson: only implemented for 1D or 3D (with a solve only in last dimension)")

   -- Solve direction is parallel (z) direction, which is always assumed to be last config space direction.
   self._zdir = self._ndim

   self._writeMatrix = xsys.pickBool(tbl.writeStiffnessMatrix, false)

   -- Set up constant dummy field.
   self.unitWeight = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   local initUnit = ProjectOnBasis {
      onGrid   = self._grid,
      basis    = self._basis,
      evaluate = function (t,xn) return 1.0 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- Set up fields for non-uniform and/or time-dependent laplacian and modifier weights
   self.laplacianWeight = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.modifierWeight = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   -- Initialize these fields to zero.
   self.laplacianWeight:clear(0.0)
   self.modifierWeight:clear(0.0)

   -- When neither laplacianWeight nor modifierWeight are set,
   -- this option effectively sets laplacian = 0 and modifier = 1.
   self._smooth = xsys.pickBool(tbl.smooth, false)
  
   local function getBcData(tbl)
      local bc = ffi.new("bcdataPar_t")
      if tbl.T == "D" then bc.type = DIRICHLET
      elseif tbl.T == "N" then bc.type = NEUMANN
      elseif tbl.T == "D_VAR" then bc.type = DIRICHLET_VARIABLE
      else assert(false, "Updater.FemParPoisson: boundary condition type must be specified by one of 'D', 'N', or 'D_VAR'")
      end
      bc.value = tbl.V
      bc.isSet = true
      return bc
   end

   -- Read in boundary conditions. 
   self._bc = ffi.new("bcdataPar_t[2]")
   self._periodic = false
   if tbl.bcLower and tbl.bcUpper then
      assert(#tbl.bcLower==1 and #tbl.bcUpper==1, "Updater.FemPerpPoisson: number of entries in bcLower/bcUpper must equal 1.")
      if tbl.bcLower[1].T=="P" and tbl.bcUpper[1].T=="P" then  
         self._periodic = true
      else
         assert(tbl.bcLower[1].T~="P" and tbl.bcUpper[1].T~="P", "Updater.FemParPoisson: both or neither 'bcLower.T' and 'bcUpper.T' have to be 'P' (periodic).")
         self._bc[0] = getBcData(tbl.bcLower[1])
         self._bc[1] = getBcData(tbl.bcUpper[1])
      end
   else
      assert(false, "Updater.FemParPoisson: must specify 'bcLower' and 'bcUpper'.")
   end

   self._hasLaplacian = false
   self._hasModifier  = false

   self._adjustSource = false

   self._nz = self._grid:numCells(self._zdir)
   self._p  = self._basis:polyOrder()
   self._dz = self._grid:dx(self._zdir)

   self._poisson   = {}
   self._first     = true
   self._makeStiff = true


   if GKYL_HAVE_MPI then
      -- Split communicators in x-y.
      local commSet   = self._grid:commSet()
      local worldComm = commSet.comm
      local nodeComm  = commSet.nodeComm
      local nodeRank  = Mpi.Comm_rank(nodeComm)
      local xyrank = 0
      if self._ndim>1 then xyrank = nodeRank%(self._grid:cuts(1)*self._grid:cuts(2)) end
      self._xycomm = Mpi.Comm_split(worldComm, xyrank, nodeRank)
   end

   self.dynVec = DataStruct.DynVector { numComponents = 1 }

   return self
end

-- For testing.
function FemParPoisson:bcType(side) return self._bc[side].type end
function FemParPoisson:bcValue(side) return self._bc[side].value end

-- Advance method.
function FemParPoisson:_advance(tCurr, inFld, outFld) 
   local grid  = self._grid
   local basis = self._basis

   local src = assert(inFld[1], "FemParPoisson.advance: must specify an input field")
   local sol = assert(outFld[1], "FemParPoisson.advance: must specify an output field")

   local ndim = self._ndim

   -- Create region that is effectively 1d and global in z directions.
   local parRange       = src:globalRange()
   local localRange     = src:localRange()
   local local_xy_lower = {1, 1}
   local local_xy_upper = {1, 1}
   for d = 1, self._ndim-1 do
      parRange:shorten(d)
      local_xy_lower[d] = localRange:lower(d)
      local_xy_upper[d] = localRange:upper(d)
   end

   -- create indexers and pointers for src and sol
   if self._first then 
      self.srcIndexer             = src:indexer() 
      self.srcPtr                 = src:get(1)
      self.solIndexer             = sol:indexer() 
      self.solPtr                 = sol:get(1)
      self.laplacianWeightIndexer = self.laplacianWeight:indexer() 
      self.laplacianWeightPtr     = self.laplacianWeight:get(1)
      self.modifierWeightIndexer  = self.modifierWeight:indexer() 
      self.modifierWeightPtr      = self.modifierWeight:get(1)
      if self._hasModifier == false and self._hasLaplacian == false then 
         if self._smooth then 
           self.modifierWeight:copy(self.unitWeight) 
           self._hasModifier = true
         else 
           self.laplacianWeight:copy(self.unitWeight)
           self._hasLaplacian = true
         end
      end
      if self._periodic and self._hasModifier == false then self._adjustSource = true end
   end

   local intSrcVol = {0.0}
   -- if all directions periodic need to adjust source so that integral is 0 
   if self._adjustSource then
      -- integrate source
      if self._first then
         self.calcInt = CartFieldIntegratedQuantCalc {
            onGrid        = grid,
            basis         = basis,
            numComponents = 1,
            quantity      = "V",
         }
      end
      self.calcInt:advance(0.0, {src}, {self.dynVec})
      _, intSrcVol = self.dynVec:lastData()
   end

   -- Loop over local x-y cells.
   for idx=local_xy_lower[1], local_xy_upper[1] do
      if self._first then
         self._poisson[idx] = {}
      end
      for idy=local_xy_lower[2], local_xy_upper[2] do
         if self._first then
            -- If first time, create poisson C object for each z cell.
            self._poisson[idx][idy] = ffiC.new_FemParPoisson(self._nz, self._ndim, self._p, 
                                               self._dz, self._periodic, 
                                               self._bc, self._writeMatrix)
         end

         -- Zero global source.
         ffiC.zeroParGlobalSrc(self._poisson[idx][idy])

         -- Create global source.
         -- globalSrc is an Eigen vector managed in C
         -- each proc allocates a full globalSrc vector
         -- loop over z cells locally to get local contributions to globalSrc
         for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
            if ndim==1 then 
               src:fill(self.srcIndexer(idz), self.srcPtr)
            elseif ndim==3 then
               src:fill(self.srcIndexer(idx, idy, idz), self.srcPtr)
            end
            ffiC.createParGlobalSrc(self._poisson[idx][idy], self.srcPtr:data(), idz-1, intSrcVol[1]/grid:gridVolume()*math.sqrt(2)^ndim)

            if self._makeStiff then
               if ndim==1 then 
                  self.laplacianWeight:fill(self.laplacianWeightIndexer(idz), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.modifierWeightIndexer(idz), self.modifierWeightPtr) 
               elseif ndim==3 then
                  self.laplacianWeight:fill(self.laplacianWeightIndexer(idx, idy, idz), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.modifierWeightIndexer(idx, idy, idz), self.modifierWeightPtr) 
               end
               ffiC.makeParGlobalStiff(self._poisson[idx][idy], self.laplacianWeightPtr:data(), self.modifierWeightPtr:data(), idz-1)
            end 
         end

         if GKYL_HAVE_MPI and Mpi.Comm_size(self._xycomm)>1 then
            -- Sum each proc's globalSrc to get final globalSrc (on each proc via allreduce).
            ffiC.allreduceParGlobalSrc(self._poisson[idx][idy], Mpi.getComm(self._xycomm))
            if self._makeStiff then
               ffiC.allgatherParGlobalStiff(self._poisson[idx][idy], Mpi.getComm(self._xycomm))
            end
         end

         if self._makeStiff then
            ffiC.finishParGlobalStiff(self._poisson[idx][idy])
         end
         -- Solve.
         ffiC.solvePar(self._poisson[idx][idy])
  
         -- Remap global nodal solution to local modal solution.
         -- Only need to loop over local proc region.
         for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
            if ndim==1 then 
               sol:fill(self.solIndexer(idz), self.solPtr)
            elseif ndim==3 then
               sol:fill(self.solIndexer(idx, idy, idz), self.solPtr)
            end
            ffiC.getSolutionPar(self._poisson[idx][idy], self.solPtr:data(), idz-1)
         end
      end
   end

   self._first = false
   -- Reset makeStiff flag to false, until stiffness matrix changes again.
   self._makeStiff = false
end

function FemParPoisson:delete()
  ffiC.delete_FemParPoisson(self._poisson)
end

function FemParPoisson:setLaplacianWeight(weight)
   self._hasLaplacian = true
   self.laplacianWeight:copy(weight)
   -- Need to remake stiffness matrix since laplacianWeight has changed.
   self._makeStiff = true
end
function FemParPoisson:setModifierWeight(weight)
   self._hasModifier = true
   self.modifierWeight:copy(weight)
   -- Need to remake stiffness matrix since modifierWeight has changed.
   self._makeStiff = true
end

function FemParPoisson:getLaplacianWeight()
   return self.laplacianWeight
end
function FemParPoisson:getModifierWeight()
   return self.modifierWeight
end

return FemParPoisson
