-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson (or rather Helmholtz) equation
--
--   div{ epsilon nabla_perp{ phi } } + beta * phi = sigma
--
-- in perpendicular directions with FEM scheme. Perpendicular directions assumed
-- to be first two configuration-space directions
--
-- We use the following terminology:
--   epsilon: Laplacian weight.
--   beta: modifier weight.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc       = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin         = require "Lib.Linalg"
local Proto       = require "Lib.Proto"
local Range       = require "Lib.Range"
local ffi         = require "ffi"
local ffiC        = ffi.C
local xsys        = require "xsys"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local FemParPoisson  = require "Updater.FemParPoisson"
local DataStruct     = require "DataStruct"
local Grid = require("Grid")
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

          int istart[8];
          int iend[8];
          int cornerstart[8];
          int cornerend[8];
      } bcdata_t;
  typedef struct FemPerpPoisson FemPerpPoisson;
  FemPerpPoisson* new_FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix, bool adjustSource);
  void delete_FemPerpPoisson(FemPerpPoisson* f);
  void makeGlobalStiff(FemPerpPoisson* f, double *laplacianWeight, double *modifierWeight, double *gxx, double *gxy, double *gyy, int idx, int idy);
  void finishGlobalStiff(FemPerpPoisson* f);
  void createGlobalSrc(FemPerpPoisson* f, double* localSrcPtr, int idx, int idy, double intSrcVol);
  void allreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm);
  void IallreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm);
  void waitForGlobalSrcReduce(FemPerpPoisson* f, MPI_Comm comm);
  void allgatherGlobalStiff(FemPerpPoisson* f, MPI_Comm comm);
  void zeroGlobalSrc(FemPerpPoisson* f);
  void solve(FemPerpPoisson* f);
  void getSolution(FemPerpPoisson* f, double* ptr, int idx, int idy);
  void getNodalSolution(FemPerpPoisson* f, double* ptr, int idx, int idy);
]]
local DIRICHLET = 0
local NEUMANN   = 1
local DIRICHLET_VARIABLE = 2
--local bcdata
--local mt = {}
--bcdata = ffi.metatype("bcdata_t", mt)

-- FEM Poisson solver updater object.
local FemPerpPoisson = Proto(UpdaterBase)

function FemPerpPoisson:init(tbl)
   FemPerpPoisson.super.init(self, tbl)

   self.xLCFS = tbl.xLCFS
   --print("In FemPerpPoisson xLCFS =", self.xLCFS)
   self._grid  = assert(tbl.onGrid, "Updater.FemPerpPoisson: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemPerpPoisson: Must specify basis functions to use using 'basis'")

   self._ndim = self._grid:ndim()

   assert(self._basis:id()=="serendipity", "Updater.FemPerpPoisson: only implemented for modal serendipity basis")
   assert(self._ndim == self._basis:ndim(), "Updater.FemPerpPoisson: dimensions of basis and grid must match")
   assert(self._basis:polyOrder()==1 or self._basis:polyOrder()==2, "Updater.FemPerpPoisson: only implemented for polyOrder = 1 or 2")
   assert(self._ndim==2 or self._ndim==3, "Updater.FemPerpPoisson: only implemented for 2D or 3D (with no solve in 3rd dimension)")

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
      onGhosts = true
   }
   initUnit:advance(0.,{},{self.unitWeight})

   -- Set up fields for non-uniform and/or time-dependent laplacian and modifier weights.
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

   self.zContinuous = xsys.pickBool(tbl.zContinuous, false)
   if self._smooth then self.zContinuous = true end

   local function getBcData(tbl)
      local bc = ffi.new("bcdata_t")
      if tbl.T == "D" then bc.type = DIRICHLET
      elseif tbl.T == "N" then bc.type = NEUMANN
      elseif tbl.T == "D_VAR" then bc.type = DIRICHLET_VARIABLE
      else assert(false, "Updater.FemPerpPoisson: boundary condition type must be specified by one of 'D', 'N', or 'D_VAR'")
      end
      bc.value = tbl.V
      bc.isSet = true
      return bc
   end

   -- Read in boundary conditions.
   self._bc               = ffi.new("bcdata_t[2][2]")
   self._allPeriodic      = true                 -- Flag to indicate all directions are periodic.
   self._isDirPeriodic    = ffi.new("bool[2]")   -- C flags to indicate if direction is periodic (0-based).
   self._isDirPeriodic[0] = false
   self._isDirPeriodic[1] = false
   if tbl.bcLower and tbl.bcUpper then
      if self._smooth or self.zContinuous then
         assert(#tbl.bcLower==self._ndim and #tbl.bcUpper==self._ndim,
                "Updater.FemPerpPoisson: number of entries in bcLower/bcUpper must equal the number of dimensions.")
      end
      for d = 1,2 do  -- Only check perpendicular BCs here.
         if tbl.bcLower[d].T=="P" and tbl.bcUpper[d].T=="P" then
            self._isDirPeriodic[d-1] = true
         else
            assert(tbl.bcLower[d].T~="P" and tbl.bcUpper[d].T~="P",
                   "Updater.FemPerpPoisson: both or neither 'bcLower.T' and 'bcUpper.T' have to be 'P' (periodic).")
            self._bc[d-1][0] = getBcData(tbl.bcLower[d])
            self._bc[d-1][1] = getBcData(tbl.bcUpper[d])
         end
         self._allPeriodic = self._allPeriodic and self._isDirPeriodic[d-1]
      end
   else
      assert(false, "Updater.FemPerpPoisson: must specify 'bcLower' and 'bcUpper'.")
   end

   self._hasLaplacian = false
   self._hasModifier  = false

   self._adjustSource = false

   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self.axisymmetric=nil
   if self._ny==1 and self.xLCFS then
     self.axisymmetric=true
     --print("yes we are axisymmetric")
   end
   self._p  = self._basis:polyOrder()
   self._dx = self._grid:dx(1)
   self._dy = self._grid:dx(2)

   if self._p==2 then assert(self._nx>1 and self._ny>1, "Updater.FemPerpPoisson: must use nx>1 and ny>1 for p=2") end

   self._poisson   = {}
   self._first     = true
   self._makeStiff = true

   -- Metric coefficient fields.
   self.gxx = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.gxy = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.gyy = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   if tbl.gxx then 
      self.gxx:copy(tbl.gxx)
   else
      self.gxx:copy(self.unitWeight)
   end
   if tbl.gxy then 
      self.gxy:copy(tbl.gxy)
   else
      self.gxy:clear(0.0)
   end
   if tbl.gyy then 
      self.gyy:copy(tbl.gyy)
   else
      self.gyy:copy(self.unitWeight)
   end

   if GKYL_HAVE_MPI then
      -- Split communicators in z.
      local commSet   = self._grid:commSet()
      local worldComm = commSet.comm
      local nodeComm  = commSet.nodeComm
      local nodeRank  = Mpi.Comm_rank(nodeComm)
      local zrank     = 0
      if self._ndim==3 then zrank = math.floor(nodeRank/self._grid:cuts(1)/self._grid:cuts(2)) end
      self._zcomm = Mpi.Comm_split(worldComm, zrank, nodeRank)
   end

   self.zDiscontToCont = nil
   if self._ndim == 3 and self.zContinuous then
      --move this up for now
      self.zDiscontToCont = FemParPoisson {
         onGrid  = self._grid,
         basis   = self._basis,
         bcLower = {tbl.bcLower[3]},
         bcUpper = {tbl.bcUpper[3]},
         smooth  = true,
      }
      if self.axisymmetric then
        --print("In FemPerpPoisson, Looking at bl.bcLower[3]: ", tbl.bcLower[3])
         self:setSplitRanges()
         self.smoother = function(tCurr,src) return self:axisymmetricSmoother(tCurr,src) end
         self.zDiscontToContSOL = FemParPoisson {
            onGrid  = self._grid,
            basis   = self._basis,
            bcLower = {{T= "N",V=0.0}},
            bcUpper = {{T= "N",V=0.0}},
            smooth  = true,
         }
         self.zDiscontToContCore = FemParPoisson {
            onGrid  = self._grid,
            basis   = self._basis,
            bcLower = {{T = "P",V=0.0}},
            bcUpper = {{T = "P",V=0.0}},
            smooth  = true,
         }
         local function createField(grid, basis, ghostCells, vComp, periodicSync)
            vComp = vComp or 1
            local fld = DataStruct.Field {
               onGrid           = grid,
               numComponents    = basis:numBasis()*vComp,
               ghost            = ghostCells,
               metaData         = {polyOrder = basis:polyOrder(),
                                   basisType = basis:id()},
               syncPeriodicDirs = periodicSync,
            }
            fld:clear(0.0)
            return fld
         end
         self.srcSOL = createField(self._grid,self._basis,{1,1})
         --self.srcCore = createField(self._grid,self._basis,{1,1})
         --self.srcTemp = createField(self._grid,self._basis,{1,1})
      else
         self.smoother = function(tCurr,src) return self:regularSmoother(tCurr,src) end
         --self.zDiscontToCont = FemParPoisson {
         --   onGrid  = self._grid,
         --   basis   = self._basis,
         --   bcLower = {tbl.bcLower[3]},
         --   bcUpper = {tbl.bcUpper[3]},
         --   smooth  = true,
         --}
      end
   end

   local localRange = self._grid:localRange()
   -- Create region that is effectively 2d and global in x-y directions.
   self.local_z_lower = 1
   self.local_z_upper = 1
   if (self._ndim == 3) then
      self.local_z_lower = localRange:lower(3)
      self.local_z_upper = localRange:upper(3)
   end

   self.dynVec = DataStruct.DynVector { numComponents = 1 }

   -- Timers.
   self.timers = {srcInt         = 0., stiffReduce    = 0.,
                  objCreate      = 0., stiffFinish    = 0.,
                  srcCreate      = 0., solve          = 0.,
                  stiffCreate    = 0., getSol         = 0.,
                  srcReduce      = 0., srcReduceWait = 0.,  assemble       = 0., 
                  completeNsolve = 0., zDiscontToCont = 0.}
   self.timerLabelsAssembly = {
                       srcInt         = "FemPerpPoisson... Assembly: Integrate source (for all-periodic BCs only)", 
                       objCreate      = "FemPerpPoisson... Assembly: Create FemPoisson object (first time only)",
                       srcCreate      = "FemPerpPoisson... Assembly: Initialize global source vector", 
                       srcReduce      = "FemPerpPoisson... Assembly: All-reduce global source vector",
                       stiffCreate    = "FemPerpPoisson... Assembly: Initialize stiffness matrix triplet list",
                       stiffReduce    = "FemPerpPoisson... Assembly: All-gather of stiffness matrix entries",
                       stiffFinish    = "FemPerpPoisson... Assembly: Assemble stiffness matrix, apply BCs, and factorize",
                       --assemble       = "FemPerpPoisson... Total assembly",
                       }
    self.timerLabelsSolve = {
                       srcReduceWait  = "FemPerpPoisson... Solve: wait for global source reduce to finish",
                       solve          = "FemPerpPoisson... Solve: solve linear system",
                       getSol         = "FemPerpPoisson... Solve: get solution",
                       zDiscontToCont = "FemPerpPoisson... Solve: do z smoothing solve",
                       --completeNsolve = "FemPerpPoisson... Total solve",
                       }

   return self
end

-- For testing.
function FemPerpPoisson:bcType(dir,side) return self._bc[dir][side].type end
function FemPerpPoisson:bcValue(dir,side) return self._bc[dir][side].value end

function FemPerpPoisson:regularSmoother(tCurr,src)
   self.zDiscontToCont:advance(tCurr, {src}, {src})
end

function FemPerpPoisson:setSplitRanges()
   local gridRange = self._grid:globalRange()
   local xLCFS=self.xLCFS
   local coordLCFS = {xLCFS-1.e-7}
   idxLCFS    = {-9}
   local xGridIngr = self._grid:childGrid({1})
   local xGrid = Grid.RectCart {
      lower = xGridIngr.lower,  periodicDirs  = xGridIngr.periodicDirs,
      upper = xGridIngr.upper,  decomposition = xGridIngr.decomposition,
      cells = xGridIngr.cells,
   }
   xGrid:findCell(coordLCFS, idxLCFS)
   self.globalCore = gridRange:shorten(1, idxLCFS[1])
   --local lv,uv = self.globalCore:lowerAsVec(), self.globalCore:upperAsVec()
   --print("Core range lv,uv", lv[1],lv[2], lv[3], uv[1],uv[2], uv[3])
   self.globalSOL = gridRange:shortenFromBelow(1, self._grid:numCells(1)-idxLCFS[1]+1)
   --local lv,uv = self.globalSOL:lowerAsVec(), self.globalSOL:upperAsVec()
   --print("SOL range lv,uv", lv[1],lv[2], lv[3], uv[1],uv[2], uv[3])
end

function FemPerpPoisson:axisymmetricSmoother(tCurr,src)
   --print("Chose the axisymmetric Smoother")
   --local function createField(grid, basis, ghostCells, vComp, periodicSync)
   --   vComp = vComp or 1
   --   local fld = DataStruct.Field {
   --      onGrid           = grid,
   --      numComponents    = basis:numBasis()*vComp,
   --      ghost            = ghostCells,
   --      metaData         = {polyOrder = basis:polyOrder(),
   --                          basisType = basis:id()},
   --      syncPeriodicDirs = periodicSync,
   --   }
   --   fld:clear(0.0)
   --   return fld
   --end
   local localRange = src:localRange()
   local localSOLRange = src:localRange():intersect(self.globalSOL)
   local localCoreRange = src:localRange():intersect(self.globalCore)
   --self:getSplitRanges() -- now called in init
   --srcSOL = createField(self._grid,self._basis,{1,1})
   --srcSOL=self.srcSOL
   --srcCore = createField(self._grid,self._basis,{1,1})
   --Have to be careful!! turns out zDiscontToCont=FemParPoisson actually modifies src so create a temporary copy
   --srcTemp = createField(self._grid,self._basis,{1,1})
   --srcTemp=self.srcTemp
   --self.srcTemp:copy(localRange, src)
   --Mpi.Barrier(self._grid:commSet().sharedComm)
   self.zDiscontToContCore:advance(tCurr, {src}, {src})
   self.zDiscontToContSOL:advance(tCurr, {src}, {self.srcSOL})
   src:copyRange(localSOLRange, self.srcSOL)
   --src:copyRange(localCoreRange, self.srcCore)
end

--function FemPerpPoisson:smootherImpl(tCurr,src)
--   return self.smoother(tCurr,src)
--end

function FemPerpPoisson:beginAssembly(tCurr, src)
   -- Assemble the right-side source vector and, if necessary, the stiffness matrix.
   local tm = Time.clock()

   local ndim, grid = self._ndim, self._grid

   --First smoothing: of density with z BCs
   if self.zDiscontToCont then   -- Make continuous in z.
      --print("doing the first overall smooth")
      local tmStart = Time.clock()
      --self.zDiscontToCont:advance(tCurr, {src}, {src})
      --call some function, function is determined at t=0
      --self.smootherImpl(tCurr,src)
      self.smoother(tCurr,src)
      self.timers.zDiscontToCont = self.timers.zDiscontToCont + Time.clock() - tmStart
   end

   local localRange = src:localRange()

   -- Create indexers and pointers for src and sol.
   if self._first then 
      self.indexer = src:indexer() 
      self.srcPtr  = src:get(1)
      self.laplacianWeightPtr = self.laplacianWeight:get(1)
      self.modifierWeightPtr  = self.modifierWeight:get(1)
      self.gxxPtr = self.gxx:get(1)
      self.gxyPtr = self.gxy:get(1)
      self.gyyPtr = self.gyy:get(1)
      if self._hasModifier == false and self._hasLaplacian == false then 
         if self._smooth then 
            self.modifierWeight:copy(self.unitWeight) 
            self._hasModifier = true
         else 
            self.laplacianWeight:copy(self.unitWeight)
            self._hasLaplacian = true
         end
      end
      if self._allPeriodic and self._hasModifier == false then self._adjustSource = true end
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

   -- Loop over local z cells.
   for idz=self.local_z_lower, self.local_z_upper do
      if self._first then 
         local tmStart = Time.clock()
         -- If first time, create poisson C object for each z cell.
         self._poisson[idz] = ffiC.new_FemPerpPoisson(self._nx, self._ny, self._ndim, self._p, 
                                                      self._dx, self._dy, self._isDirPeriodic, 
                                                      self._bc, self._writeMatrix, self._adjustSource)
         self.timers.objCreate = self.timers.objCreate + Time.clock() - tmStart
      end

      -- Zero global source.
      ffiC.zeroGlobalSrc(self._poisson[idz])

      -- Create global source.
      -- globalSrc is an Eigen vector managed in C.
      -- Each proc allocates a full globalSrc vector.
      -- Loop over x and y cells locally to get local contributions to globalSrc.
      for idx=localRange:lower(1),localRange:upper(1) do     
         for idy=localRange:lower(2),localRange:upper(2) do
            local tmStart = Time.clock()
            if ndim==2 then 
               src:fill(self.indexer(idx,idy), self.srcPtr) 
            else 
               src:fill(self.indexer(idx, idy, idz), self.srcPtr) 
            end
            ffiC.createGlobalSrc(self._poisson[idz], self.srcPtr:data(), idx-1, idy-1, intSrcVol[1]/grid:gridVolume()*math.sqrt(2)^self._ndim)
            self.timers.srcCreate = self.timers.srcCreate + Time.clock() - tmStart

            if self._makeStiff then
               local tmStiffStart = Time.clock()
               if ndim==2 then 
                  self.laplacianWeight:fill(self.indexer(idx,idy), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.indexer(idx,idy), self.modifierWeightPtr) 
                  self.gxx:fill(self.indexer(idx,idy), self.gxxPtr) 
                  self.gxy:fill(self.indexer(idx,idy), self.gxyPtr) 
                  self.gyy:fill(self.indexer(idx,idy), self.gyyPtr) 
               else 
                  self.laplacianWeight:fill(self.indexer(idx, idy, idz), self.laplacianWeightPtr) 
                  self.modifierWeight:fill(self.indexer(idx, idy, idz), self.modifierWeightPtr) 
                  self.gxx:fill(self.indexer(idx,idy,idz), self.gxxPtr) 
                  self.gxy:fill(self.indexer(idx,idy,idz), self.gxyPtr) 
                  self.gyy:fill(self.indexer(idx,idy,idz), self.gyyPtr) 
               end
               ffiC.makeGlobalStiff(self._poisson[idz], self.laplacianWeightPtr:data(), self.modifierWeightPtr:data(),
                                    self.gxxPtr:data(), self.gxyPtr:data(), self.gyyPtr:data(), idx-1, idy-1)
               self.timers.stiffCreate = self.timers.stiffCreate + Time.clock() - tmStiffStart
            end 
         end
      end

      if GKYL_HAVE_MPI and Mpi.Comm_size(self._zcomm)>1 then
         local tmStart = Time.clock()
         -- Sum each proc's globalSrc to get final globalSrc (on each proc via allreduce).
         ffiC.IallreduceGlobalSrc(self._poisson[idz], Mpi.getComm(self._zcomm))
         self.timers.srcReduce = self.timers.srcReduce + Time.clock() - tmStart
         if self._makeStiff then
            local tmStiffStart = Time.clock()
            ffiC.allgatherGlobalStiff(self._poisson[idz], Mpi.getComm(self._zcomm))
            self.timers.stiffReduce = self.timers.stiffReduce + Time.clock() - tmStiffStart
         end
      end

      if self._makeStiff then
         local tmStiffStart = Time.clock()
         ffiC.finishGlobalStiff(self._poisson[idz])
         self.timers.stiffFinish = self.timers.stiffFinish + Time.clock() - tmStiffStart
      end
   end

   self.timers.assemble = self.timers.assemble + Time.clock() - tm
end

function FemPerpPoisson:completeAssemblyAndSolve(tCurr, sol)
   -- Wait for the reduction of the source vector and (if needed) the stiffness
   -- matrix, and perform the linear solve.
   local tm = Time.clock()

   local ndim       = self._ndim
   local localRange = sol:localRange()

   if self._first then self.solPtr = sol:get(1) end

   -- Loop over local z cells.
   for idz=self.local_z_lower, self.local_z_upper do
      if GKYL_HAVE_MPI and Mpi.Comm_size(self._zcomm)>1 then
         local tmStart = Time.clock()
         -- Wait for the sum of each proc's globalSrc to get final globalSrc.
         ffiC.waitForGlobalSrcReduce(self._poisson[idz], Mpi.getComm(self._zcomm))
         self.timers.srcReduceWait = self.timers.srcReduceWait + Time.clock() - tmStart
      end

      -- Solve. 
      local tmSolStart = Time.clock()
      ffiC.solve(self._poisson[idz])
      self.timers.solve = self.timers.solve + Time.clock() - tmSolStart

      -- Remap global nodal solution to local modal solution.
      -- Only need to loop over local proc region.
      local tmGetStart = Time.clock()
      for idx=localRange:lower(1),localRange:upper(1) do     
         for idy=localRange:lower(2),localRange:upper(2) do
            if ndim==2 then 
               sol:fill(self.indexer(idx,idy), self.solPtr) 
            else 
               sol:fill(self.indexer(idx, idy, idz), self.solPtr) 
            end
            ffiC.getSolution(self._poisson[idz], self.solPtr:data(), idx-1, idy-1)
         end
      end
      self.timers.getSol = self.timers.getSol + Time.clock() - tmGetStart
   end 

   --Second Smoothing: smoothing of phi with Z BCs
   --if self.zDiscontToCont then 
   --   print("doing the second overall smooth")
   --   local tmStart = Time.clock()
   --   self.zDiscontToCont:advance(tCurr, {sol}, {sol}) 
   --   self.timers.zDiscontToCont = self.timers.zDiscontToCont + Time.clock() - tmStart
   --end

   self._first = false
   -- Reset makeStiff flag to false, until stiffness matrix changes again.
   self._makeStiff = false

   self.timers.completeNsolve = self.timers.completeNsolve + Time.clock() - tm
end

function FemPerpPoisson:assemble(tCurr, inFld, outFld)
   -- Begin assembling the source vector and, if needed, the stiffness matrix.
   local src = assert(inFld[1], "FemPerpPoisson.advance: Must specify an input field")

   self:beginAssembly(tCurr, src)
end

function FemPerpPoisson:solve(tCurr, inFld, outFld)
   -- Assuming the linear problem has been assemble (or the assembly has been initiated
   -- via non-blocking MPI), :solve waits for the assembly to finish and solves the problem.
   local sol = assert(outFld[1], "FemPerpPoisson.advance: Must specify an output field")

   -- Compute and obtain the local solution.
   self:completeAssemblyAndSolve(tCurr, sol)
end

function FemPerpPoisson:_advance(tCurr, inFld, outFld) 
   -- Advance method. Assembles the linear problem and solves it.
   local ndim  = self._ndim

   local src = assert(inFld[1], "FemPerpPoisson.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemPerpPoisson.advance: Must specify an output field")

   -- Assemble the right-side source vector and, if needed, the stiffness matrix.
   self:beginAssembly(tCurr, src)

   -- Compute and obtain the local solution.
   self:completeAssemblyAndSolve(tCurr, sol)
end

function FemPerpPoisson:setLaplacianWeight(weight)
   self._hasLaplacian = true
   self.laplacianWeight:copy(weight)
   -- Need to remake stiffness matrix since laplacianWeight has changed.
   self._makeStiff = true
end

function FemPerpPoisson:setModifierWeight(weight)
   self._hasModifier = true
   self.modifierWeight:copy(weight)
   -- Need to remake stiffness matrix since modifierWeight has changed.
   self._makeStiff = true
end

function FemPerpPoisson:getLaplacianWeight()
   return self.laplacianWeight
end

function FemPerpPoisson:getModifierWeight()
   return self.modifierWeight
end

function FemPerpPoisson:delete()
   for idz=self.local_z_lower, self.local_z_upper do
      ffiC.delete_FemPerpPoisson(self._poisson[idz])
   end
end

function FemPerpPoisson:printDevDiagnostics()
  -- Print performance/numerical diagnostics.
  local log = Logger{logToFile = true}
  log("\n")
  local calcTimeTot = math.max(self.totalTime,self.timers.assemble+self.timers.completeNsolve)
  for nm, str in pairs(self.timerLabelsAssembly) do
     v = self.timers[nm]
     log(string.format("%-99s %12.5f sec  (%6.3f%%) \n", str .. " took", v, 100*v/calcTimeTot))
  end
  for nm, str in pairs(self.timerLabelsSolve) do
     v = self.timers[nm]
     log(string.format("%-99s %12.5f sec  (%6.3f%%) \n", str .. " took", v, 100*v/calcTimeTot))
  end
  log("\n")
  v = self.timers.assemble
  log(string.format("%-99s %12.5f sec  (%6.3f%%) \n", "FemPerpPoisson... Total assembly took", v, 100*v/calcTimeTot))
  v = self.timers.completeNsolve
  log(string.format("%-99s %12.5f sec  (%6.3f%%) \n", "FemPerpPoisson... Total solve took", v, 100*v/calcTimeTot))

  log("\n")
  v = calcTimeTot
  log(string.format("%-99s %12.5f sec  (%6.3f%%) \n", "FemPerpPoisson... Total for advance", v, 100*v/calcTimeTot))
end

return FemPerpPoisson
