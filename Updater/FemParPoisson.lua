-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson (or rather Helmholtz) equation 
--
--   partial_z{ epsilon partial_z{ phi } } + beta * phi = sigma
--
-- in parallel direction with FEM scheme. Parallel direction assumed to be last
-- configuration-space direction (z).
-- 
-- We use the following terminology:
--   epsilon: Laplacian weight.
--   beta: modifier weight.
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
local Time           = require "Lib.Time"
local Logger         = require "Lib.Logger"
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
  void IallreduceParGlobalSrc(FemParPoisson* f, MPI_Comm comm);
  void waitForParGlobalSrcReduce(FemParPoisson* f, MPI_Comm comm);
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

   if self._ndim>1 then
      assert(self._basis:id()=="serendipity", "Updater.FemParPoisson: only implemented for modal serendipity basis")
   end
   assert(self._basis:polyOrder()==1 or self._basis:polyOrder()==2, "Updater.FemParPoisson: only implemented for polyOrder = 1 or 2")
   assert(self._ndim == self._basis:ndim(), "Updater.FemParPoisson: Dimensions of basis and grid must match")
   assert(self._ndim==1 or self._ndim==3, "Updater.FemParPoisson: only implemented for 1D or 3D (with a solve only in last dimension)")

   -- Solve direction is parallel (z) direction, which is always assumed to be last config space direction.
   self._zdir = self._ndim

   self._writeMatrix = xsys.pickBool(tbl.writeStiffnessMatrix, false)

   -- Set up constant dummy field.
   self.unitWeight = DataStruct.Field {
      onGrid        = self._grid,              ghost     = {1, 1},
      numComponents = self._basis:numBasis(),  useDevice = false,
   }
   local initUnit = ProjectOnBasis {
      onGrid = self._grid,   evaluate = function (t,xn) return 1.0 end,
      basis  = self._basis,  onGhosts = true,
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- Set up fields for non-uniform and/or time-dependent laplacian and modifier weights
   self.laplacianWeight = DataStruct.Field {
      onGrid        = self._grid,              ghost     = {1, 1},
      numComponents = self._basis:numBasis(),  useDevice = false,
   }
   self.modifierWeight = DataStruct.Field {
      onGrid        = self._grid,              ghost     = {1, 1}, 
      numComponents = self._basis:numBasis(),  useDevice = false,
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
      assert(#tbl.bcLower==1 and #tbl.bcUpper==1, "Updater.FemParPoisson: number of entries in bcLower/bcUpper must equal 1.")
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

   -- Create region that is effectively 1D and global in z directions.
   self.local_xy_lower = {1, 1}
   self.local_xy_upper = {1, 1}
   local localRange    = self._grid:localRange()
   for d = 1, self._ndim-1 do
      self.local_xy_lower[d] = localRange:lower(d)
      self.local_xy_upper[d] = localRange:upper(d)
   end


   -- Timers.
   self.timers = {srcInt         = 0., stiffReduce = 0.,
                  objCreate      = 0., stiffFinish = 0.,
                  srcCreate      = 0., solve       = 0.,
                  stiffCreate    = 0., getSol      = 0.,
                  srcReduce      = 0., assemble    = 0.,
                  completeNsolve = 0.}

   return self
end

-- For testing.
function FemParPoisson:bcType(side) return self._bc[side].type end
function FemParPoisson:bcValue(side) return self._bc[side].value end

function FemParPoisson:beginAssembly(tCurr, src)
   -- Assemble the right-side source vector and, if necessary, the stiffness matrix.
   local tm = Time.clock()

   local ndim, grid = self._ndim, self._grid

   if self._first then 
      -- Create indexers and pointers, and set some flags.
      self.indexer = src:indexer() 
      self.srcPtr  = src:get(1)
      self.laplacianWeightPtr = self.laplacianWeight:get(1)
      self.modifierWeightPtr  = self.modifierWeight:get(1)
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
   -- If all directions periodic need to adjust source so that integral is 0.
   if self._adjustSource then
      local tmStart = Time.clock()
      -- Integrate source.
      if self._first then
         self.calcInt = CartFieldIntegratedQuantCalc {
            onGrid        = grid,
            basis         = self._basis,
            numComponents = 1,
            quantity      = "V",
         }
      end
      self.calcInt:advance(0.0, {src}, {self.dynVec})
      _, intSrcVol = self.dynVec:lastData()
      self.timers.srcInt = self.timers.srcInt + Time.clock() - tmStart 
   end

   local localRange = src:localRange()

   -- Loop over local x-y cells.
   for idx=self.local_xy_lower[1], self.local_xy_upper[1] do
      if self._first then self._poisson[idx] = {} end
      for idy=self.local_xy_lower[2], self.local_xy_upper[2] do
         if self._first then
            -- If first time, create poisson C object for each z cell.
            local tmStart = Time.clock()
            self._poisson[idx][idy] = ffiC.new_FemParPoisson(self._nz, self._ndim, self._p, 
                                                             self._dz, self._periodic, 
                                                             self._bc, self._writeMatrix)
            self.timers.objCreate = self.timers.objCreate + Time.clock() - tmStart 
         end

         -- Zero global source.
         ffiC.zeroParGlobalSrc(self._poisson[idx][idy])

         -- Create global source.
         -- globalSrc is an Eigen vector managed in C.
         -- Each proc allocates a full globalSrc vector.
         -- Loop over z cells locally to get local contributions to globalSrc.
         for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
            local tmStart = Time.clock()
            if ndim==1 then 
               src:fill(self.indexer(idz), self.srcPtr)
            elseif ndim==3 then
               src:fill(self.indexer(idx, idy, idz), self.srcPtr)
            end
            ffiC.createParGlobalSrc(self._poisson[idx][idy], self.srcPtr:data(), idz-1, intSrcVol[1]/grid:gridVolume()*math.sqrt(2)^ndim)
            self.timers.srcCreate = self.timers.srcCreate + Time.clock() - tmStart 

            if self._makeStiff then
               local tmStiffStart = Time.clock()
               if ndim==1 then 
                  self.laplacianWeight:fill(self.indexer(idz), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.indexer(idz), self.modifierWeightPtr) 
               elseif ndim==3 then
                  self.laplacianWeight:fill(self.indexer(idx, idy, idz), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.indexer(idx, idy, idz), self.modifierWeightPtr) 
               end
               ffiC.makeParGlobalStiff(self._poisson[idx][idy], self.laplacianWeightPtr:data(), self.modifierWeightPtr:data(), idz-1)
               self.timers.stiffCreate = self.timers.stiffCreate + Time.clock() - tmStiffStart 
            end 
         end

         if GKYL_HAVE_MPI and Mpi.Comm_size(self._xycomm)>1 then
            local tmStart = Time.clock()
            -- Sum each proc's globalSrc to get final globalSrc (on each proc via allreduce).
            ffiC.IallreduceParGlobalSrc(self._poisson[idx][idy], Mpi.getComm(self._xycomm))
            self.timers.srcReduce = self.timers.srcReduce + Time.clock() - tmStart 
            if self._makeStiff then
               local tmStiffStart = Time.clock()
               ffiC.allgatherParGlobalStiff(self._poisson[idx][idy], Mpi.getComm(self._xycomm))
               self.timers.stiffReduce = self.timers.stiffReduce + Time.clock() - tmStiffStart 
            end
         end

         if self._makeStiff then
            local tmStiffStart = Time.clock()
            ffiC.finishParGlobalStiff(self._poisson[idx][idy])
            self.timers.stiffFinish = self.timers.stiffFinish + Time.clock() - tmStiffStart 
         end
      end
   end

   self.timers.assemble = self.timers.assemble + Time.clock() - tm
end

function FemParPoisson:completeAssemblyAndSolve(tCurr, sol)
   -- Wait for the reduction of the source vector and (if needed) the stiffness
   -- matrix, and perform the linear solve.
   local tm = Time.clock()

   local ndim       = self._ndim
   local localRange = sol:localRange()

   if self._first then self.solPtr = sol:get(1) end

   -- Perform the linear solve and get the local solution.
   for idx=self.local_xy_lower[1], self.local_xy_upper[1] do
      for idy=self.local_xy_lower[2], self.local_xy_upper[2] do
         if GKYL_HAVE_MPI and Mpi.Comm_size(self._xycomm)>1 then
            local tmStart = Time.clock()
            -- Wait for the sum of each proc's globalSrc to get final globalSrc.
            ffiC.waitForParGlobalSrcReduce(self._poisson[idx][idy], Mpi.getComm(self._xycomm))
            self.timers.srcReduce = self.timers.srcReduce + Time.clock() - tmStart 
         end

         -- Solve.
         local tmSolStart = Time.clock()
         ffiC.solvePar(self._poisson[idx][idy])
         self.timers.solve = self.timers.solve + Time.clock() - tmSolStart 
  
         -- Remap global nodal solution to local modal solution.
         -- Only need to loop over local proc region.
         local tmGetStart = Time.clock()
         for idz = localRange:lower(self._zdir), localRange:upper(self._zdir) do
            if ndim==1 then 
               sol:fill(self.indexer(idz), self.solPtr)
            elseif ndim==3 then
               sol:fill(self.indexer(idx, idy, idz), self.solPtr)
            end
            ffiC.getSolutionPar(self._poisson[idx][idy], self.solPtr:data(), idz-1)
         end
         self.timers.getSol = self.timers.getSol + Time.clock() - tmGetStart 
      end
   end

   self._first = false
   -- Reset makeStiff flag to false, until stiffness matrix changes again.
   self._makeStiff = false

   self.timers.completeNsolve = self.timers.completeNsolve + Time.clock() - tm
end

function FemParPoisson:assemble(tCurr, inFld, outFld)
   -- Begin assembling the source vector and, if needed, the stiffness matrix.
   local src = assert(inFld[1], "FemParPoisson.advance: Must specify an input field")

   self:beginAssembly(tCurr, src)
end

function FemParPoisson:solve(tCurr, inFld, outFld)
   -- Assuming the linear problem has been assemble (or the assembly has been initiated
   -- via non-blocking MPI), :solve waits for the assembly to finish and solves the problem.
   local sol = assert(outFld[1], "FemParPoisson.advance: Must specify an output field")

   -- Compute and obtain the local solution.
   self:completeAssemblyAndSolve(tCurr, sol)
end

function FemParPoisson:_advance(tCurr, inFld, outFld) 
   -- Advance method. Assembles the linear problem and solves it.
   local src = assert(inFld[1], "FemParPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemParPoisson.advance: Must specify an output field")

   -- Assemble the right-side source vector and, if needed, the stiffness matrix.
   self:beginAssembly(tCurr, src)

   -- Compute and obtain the local solution.
   self:completeAssemblyAndSolve(tCurr, sol)
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

function FemParPoisson:delete()
  for idx=self.local_xy_lower[1], self.local_xy_upper[1] do
     for idy=self.local_xy_lower[2], self.local_xy_upper[2] do
        ffiC.delete_FemParPoisson(self._poisson[idx][idy])
     end
  end
end

function FemParPoisson:printDevDiagnostics()
  -- Print performance/numerical diagnostics.
  local log = Logger{logToFile = true}
  log("\n")
  local calcTimeTot = math.max(self.totalTime,self.timers.assemble+self.timers.completeNsolve)
  for nm, v in pairs(self.timers) do
     log(string.format(
        "FemParPoisson's "..nm.." took			%9.5f sec   (%6.3f%%)\n", v, 100*v/calcTimeTot))
  end
  log(string.format(
     "FemParPoisson's advance took			%9.5f sec   (%6.3f%%)\n",
     calcTimeTot, 100*calcTimeTot/calcTimeTot))
end

return FemParPoisson
