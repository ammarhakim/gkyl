-- Gkyl ------------------------------------------------------------------------
--
-- Constant diffusion equation on a rectangular mesh.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local Lin    = require "Lib.Linalg"
local Proto  = require "Lib.Proto"
local ffi    = require "ffi"
local xsys   = require "xsys"
local EqBase = require "Eq.EqBase"
local lume   = require "Lib.lume"
local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"
local LinearDecomp          = require "Lib.LinearDecomp"

-- ConstDiffusion equation on a rectangular mesh.
local ConstDiffusion = Proto(EqBase)

local solverTypes = { "DG", "FEM" }
local function isSolverTypeGood(typeIn)
   idxFound = lume.find(solverTypes, typeIn)
   if idxFound then return true end
   return false
end

-- Boundary condition ID numbers.
local BVP_BC_PERIODIC  = 0
local BVP_BC_DIRICHLET = 1
local BVP_BC_NEUMANN   = 2
local BVP_BC_ROBIN     = 3

-- ctor
function ConstDiffusion:init(tbl)

   self._basis = assert(tbl.basis,
      "Eq.constDiffusion: Must specify basis functions to use using 'basis'")

   solverType = tbl.solverType and tbl.solverType or "DG"
   assert(isSolverTypeGood(solverType), "Updater.ConstDiffusion: solverType must be one of 'DG' or 'FEM'.")
   if solverType == "DG" then
      self.isDG  = true
      self.isFEM = false
   elseif solverType == "FEM" then
      self.isDG  = false
      self.isFEM = true
   end

   local basisID, polyOrder = self._basis:id(), self._basis:polyOrder()
   self.dim = self._basis:ndim()

   -- Read the directions in which to apply diffusion. 
   local diffDirsIn = tbl.diffusiveDirs
   if diffDirsIn then
      assert(#diffDirsIn<=self.dim, "Eq.constDiffusion: 'diffusiveDirs' cannot have more entries than the simulation's self.dimensions.")
      if self.isFEM then
         assert(#diffDirsIn==self.dim, "Eq.constDiffusion: 'diffusiveDirs' must have the same number of entries as the simulation's self.dimensions if solverType=FEM.")
      end
      self.diffDirs = diffDirsIn
   else
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, self.dim do self.diffDirs[d] = d end
   end
   
   -- Read diffusion coefficient (or vector).
   local nuIn     = assert(tbl.coefficient,
                           "Eq.constDiffusion: must specify diffusion coefficient (or vector) using 'coefficient' ")
   local nuInType = type(nuIn)
   self._nu       = Lin.Vec(self.dim)
   for d = 1, self.dim do self._nu[d] = 0.0 end
   if (nuInType == "number") then
      -- Set the diffusion coefficient to the same amplitude in all directions.
      for d = 1, self.dim do self._nu[d] = nuIn end
   elseif (nuInType == "table") then
      if diffDirsIn then
         assert(#nuIn==#diffDirsIn, "Eq.constDiffusion: 'coefficient' table must have the same number of entries as 'diffusiveDirs'.")
      else
         assert(#nuIn==self.dim, "Eq.constDiffusion: 'coefficient' table must have the same number of entries as the simulation's self.dimensions if 'diffusiveDirs' is not given.")
      end
      for d = 1, #self.diffDirs do self._nu[self.diffDirs[d]] = nuIn[d] end
   else
      assert(false, "Eq.constDiffusion: 'coefficient' must be a number or a table.")
   end

   local diffOrder
   if tbl.order then
      diffOrder = tbl.order
      assert(not (diffOrder > 4 and polyOrder < 2), "Eq.constDiffusion: grad^6 requires polyOrder > 1.")
   else
      diffOrder = 2
   end

   local applyPositivity = xsys.pickBool(tbl.positivity, false)   -- Positivity preserving option.

   -- Store pointers to C kernels implementing volume and surface terms.
   self._volTerm           = ConstDiffusionModDecl.selectVol(basisID, self.dim, polyOrder, self.diffDirs, diffOrder)
   self._surfTerms         = ConstDiffusionModDecl.selectSurf(basisID, self.dim, polyOrder, self.diffDirs, diffOrder, applyPositivity)
   self._boundarySurfTerms = ConstDiffusionModDecl.selectBoundarySurf(basisID, self.dim, polyOrder, self.diffDirs, diffOrder, applyPositivity)

   -- Select the right volume and surface terms for FEM/DG.
   self.volTermFunc  = function(w, dx, idx, f, out) return ConstDiffusion["volTerm" .. solverType](self, w, dx, idx, f, out) end
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return ConstDiffusion["surfTerm" .. solverType](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end
   self.boundarySurfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return ConstDiffusion["boundarySurfTerm" .. solverType](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end

   if self.isFEM then
      self.grid = assert(tbl.onGrid, "Eq.ConstDiffusion (solverType=FEM): Must provide grid object using 'onGrid'.")

      local bcLower = assert(
         tbl.bcLower, "Updater.ConstDiffusion: Must provide lower-boundary BCs along each direction in 'bcLower'.")
      local bcUpper = assert(
         tbl.bcUpper, "Updater.ConstDiffusion: Must provide upper-boundary BCs along each direction in 'bcUpper'.")
      assert(#bcLower == self.dim, "Updater.ConstDiffusion: number of lower BCs must equal number of dimensions.")
      assert(#bcUpper == self.dim, "Updater.ConstDiffusion: number of lower BCs must equal number of dimensions.")
      local bcTypes  = {}
      local bcValues = {}
      for i = 1, self.dim do
         bcTypes[i]  = {bcID(bcLower[i]["T"]), bcID(bcUpper[i]["T"])}
         bcValues[i] = {bcLower[i]["V"], bcUpper[i]["V"]}
         for j = 1,2 do
            -- Ensure that for Robin BCs three values are given.
            if ( (bcTypes[i][j] == BVP_BC_ROBIN) and ((type(bcValues[i][j]) == "number") or
                 ((type(bcValues[i][j]) == "table") and (#bcValues[i][j] ~= 3))) ) then
               assert(false, "Updater.ConstDiffusion: for Robin BC please pass a table with 3 numbers in order to impose: bcValue1*f+bcValue2*df/dx=bcValue3.")
            end
            -- Turn BC value into a table if not Robin BC so they all are treated as tables below.
            if (type(bcValues[i][j]) == "number") then
               bcValues[i][j] = {bcValues[i][j]}
            end
         end
      end
      -- Establish periodic directions.
      self.periodicDirs   = {}
      local isDirPeriodic = {}
      for d = 1, self.dim do
         if ((bcTypes[d][1] == 0) and (bcTypes[d][2] == 0)) then
            lume.push(self.periodicDirs,d)
            isDirPeriodic[d] = true
            bcValues[d]      = {{0.0}, {0.0}}   -- Not used, but a nil could cause problems.
         elseif ( ((bcTypes[d][1] == 0) and (bcTypes[d][2] ~= 0)) or
                  ((bcTypes[d][1] ~= 0) and (bcTypes[d][2] == 0)) ) then
            assert(false, "Updater.ConstDiffusion: lower an upper BCs must both be T=\"P\" if periodic is desired.")
         else
            isDirPeriodic[d] = false
         end
      end
      self.isPeriodicDomain = lume.all(isDirPeriodic)
      self.aPeriodicDir     = lume.any(isDirPeriodic)

      -- Translate bcValues to a vector from which we can pass a pointer.
      -- This vector has 3 values for each boundary in order to support a Robin BC like:
      --     bcValue1*f+bcValue2*df/dx=bcValue3.
      self.bcValue = Lin.Vec(self.dim*2*3)
      for d = 1,self.dim do
         off = (d-1)*6
         for i = 1, 6 do
            self.bcValue[off+i] = 1.0
         end
         self.bcValue[off + 3 - (bcTypes[d][1] % 3)] = 0.0
         self.bcValue[off + 3]                       = bcValues[d][1][math.floor(1./3.)*2+1]
         if bcTypes[d][1] == BVP_BC_ROBIN then
            -- Robin BCs. First two values multiply boundary value and derivative, respectively.
            self.bcValue[off + 1] = bcValues[d][1][1]
            self.bcValue[off + 2] = bcValues[d][1][2]
         end
   
         self.bcValue[off + 6 - (bcTypes[d][2] % 3)] = 0.0
         self.bcValue[off + 6]                       = bcValues[d][2][math.floor(1./3.)*2+1]
         if bcTypes[d][2] == BVP_BC_ROBIN then
            -- Robin BCs. First two values multiply boundary value and derivative, respectively.
            self.bcValue[off + 4] = bcValues[d][2][1]
            self.bcValue[off + 5] = bcValues[d][2][2]
         end
      end

      self.cellsN = {}
      for d = 1, self.dim do self.cellsN[d]=self.grid:numCells(d) end

      self.zeros, self.ones  = {}, {}
      self.twos, self.threes = {}, {}
      self.mOnes = {}
      for d = 1, self.dim do
         self.zeros[d], self.ones[d]  = 0, 1
         self.twos[d], self.threes[d] = 2, 3
         self.mOnes[d] = -1
      end

      -- Select kernels for relaxation and computing the residual.
      self._diffusionFEM = ConstDiffusionModDecl.selectFEMdiff(basisID, self.dim, polyOrder, bcTypes)

      self.phiStencilWidth = 3
      -- 'Filled' stencils for FEM (see opStencilIndices).
      self.phiStencilSize = self.phiStencilWidth^self.dim
      -- FEM uses the right-side source in neighboring cells as well.
      self.phiStencilType = {2, self.threes, self.zeros}

      local DoublePtrVec = Lin.new_vec_ct(ffi.typeof("double*"))
      -- List of cell indices pointed to by the stencils.
      self.phiStencilIdx = {}
      for i = 1, self.phiStencilSize do self.phiStencilIdx[i] = Lin.IntVec(self.dim) end
      -- List of pointers to the data in cells pointed to by the stencil.
      self.phiStencil = DoublePtrVec(self.phiStencilSize)
      self.dxStencil  = DoublePtrVec(self.phiStencilSize)
   
      self.dxBuf = Lin.Vec(self.dim)       -- Buffer to store cell lengths.

      -- Dimensions remaining when a dimension is removed.
      self.dimRemain = {}
      for d1 = 1, self.dim do
         self.dimRemain[d1] = {}
         for d2 = 1, self.dim do self.dimRemain[d1][d2] = d2 end
         table.remove(self.dimRemain[d1],d1)
      end
   end
end

-- Methods.
function ConstDiffusion:numEquations() return 1 end
function ConstDiffusion:numWaves() return 1 end
function ConstDiffusion:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- Flux in direction dir.
function ConstDiffusion:flux(dir, qIn, fOut)
   assert(false, "ConstDiffusion:flux: NYI!")
end

-- Riemann problem for ConstDiffusion equation.
function ConstDiffusion:rp(dir, delta, ql, qr, waves, s)
   assert(false, "ConstDiffusion:rp: NYI!")
end

-- Compute q-fluctuations.
function ConstDiffusion:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "ConstDiffusion:qFluctuations: NYI!")
end

-- Maximum wave speed.
function ConstDiffusion:maxSpeed(dir, w, dx, q)
   assert(false, "ConstDiffusion:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term.
function ConstDiffusion:volTerm(w, dx, idx, f, out)
   return self.volTermFunc(w, dx, idx, f, out)
end
function ConstDiffusion:volTermDG(w, dx, idx, q, out)
   return self._volTerm(w:data(), dx:data(), self._nu:data(), q:data(), out:data())
end
function ConstDiffusion:volTermFEM(w, dx, idx, q, out)
   -- q: FEM potential in cell idx.
   self:diffusionTermFEM(idx, q, out)
   return 0
end

-- Surface integral term.
function ConstDiffusion:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function ConstDiffusion:surfTermDG(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end
function ConstDiffusion:surfTermFEM(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   return 0
end

-- Contribution from surface integral term at the boundaries.
function ConstDiffusion:boundarySurfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.boundarySurfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function ConstDiffusion:boundarySurfTermDG(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end
function ConstDiffusion:boundarySurfTermFEM(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end

function ConstDiffusion:setAuxFields(auxFields)
   if self.isFEM then
      self.phiFld = auxFields[1]
   end
end

--
-- Functions needed for FEM solve below.
-- Some of these are modeled after MGpoisson.
--

function ConstDiffusion:opStencilIndices(idxIn, stencilType, stencilIdx)
   -----------------
   -- Given the index of the current cell (idxIn), return a table
   -- of indices of the cells in a stencil of type 'stencilType'
   -- Stencil types are given by [t,w,l]:
   --   t: type of stencil
   --        0 : 'cross'.
   --        1 : partially filled cross (include nearest corner cells).
   --        2 : filled (include all corner cells).
   --   w: width of the stencil (along each dimension).
   --   l: location in grid along each dimension. For 1D
   --        [-1] = lower boundary cell.
   --        [ 0] = inner cell.
   --        [ 1] = upper boundary cell.
   --      Note that depending on how the boundary kernels are codded
   --      it may be ok to pass an inner-cell stencil at the boundary.
   --
   -- Examples:
   -- a) stencilType=[0,[5,5],[0,0]] corresponds to
   --                 9
   --                 7
   --         4   2   1   3   5
   --                 6
   --                 8
   -- b) stencilType=[1,[3,3],[0,0]] or [2,[3,3],[0,0]] correspond to
   --             7   5   9
   --             2   1   3
   --             6   4   8

   -- First copy all indicies since (in higher dimensions)
   -- most of the stay the same for each cell.
   for i = 1, #stencilIdx do
      idxIn:copyInto(stencilIdx[i])
   end

   if stencilType[1] == 0 then
      local sI = 1
      for d = 1, self.dim do
         for pm = 1,stencilType[2][d]-1 do
            sI = sI + 1
            stencilIdx[sI][d] = idxIn[d]+((-1)^(pm % 2))*((stencilType[2][d]-1)/2)
         end
      end
   elseif stencilType[1]==2 then
      local sI = 1
      for d = 1, self.dim do
         local prevDimCells = sI
         for pDC = 1, prevDimCells do
            for pm = 1-stencilType[3][d],stencilType[2][d]-1 do
               sI = sI + 1
               for _, dr in ipairs(self.dimRemain[d]) do stencilIdx[sI][dr] = stencilIdx[pDC][dr] end
               stencilIdx[sI][d] = stencilIdx[pDC][d]
                                  +((-1)^(pm % 2))*((stencilType[2][d]-1)/(1+stencilType[2][d]-stencilType[1]))
            end
         end
      end
   end
end

function ConstDiffusion:diffusionTermFEM(idx, phiItr, outItr)
   -- Compute the diffusion term, nabla^{2m}(phi), using FEM.
   local grid = self.grid

   grid:setIndex(idx)

   -- Get with indices of cells used by stencil. Store them in self.phiStencilIdx.
   self:opStencilIndices(idx, self.phiStencilType, self.phiStencilIdx)

   if (self.isPeriodicDomain and (self.dim > 1)) then
      -- Kernels need the corner ghost cells. Flip cell index to opposite corner.
      local checkFor, setTo = {{},{}}, {{},{}}
      for d = 1, self.dim do
         checkFor[1][d] = 0
         setTo[1][d]    = self.cellsN[d]
         checkFor[2][d] = self.cellsN[d]+1
         setTo[2][d]    = 1
      end
      for nI = 1, #self.phiStencilIdx do
         local changeIdxInDir = {}
         for d = 1, self.dim do
            if (self.phiStencilIdx[nI][d] == checkFor[1][d]) then  -- Lower ghost cell in this direction.
               for _, dr in ipairs(self.dimRemain[d]) do
                  if (self.phiStencilIdx[nI][dr] == checkFor[1][dr]) or
                     (self.phiStencilIdx[nI][dr] == checkFor[2][dr]) then  -- Lower/upper ghost cell in another direction too.
                     table.insert(changeIdxInDir, {d, setTo[1][d]})
                     break
                  end
               end
            elseif (self.phiStencilIdx[nI][d] == checkFor[2][d]) then  -- Lower ghost cell in this direction.
               for _, dr in ipairs(self.dimRemain[d]) do
                  if (self.phiStencilIdx[nI][dr] == checkFor[1][dr]) or
                     (self.phiStencilIdx[nI][dr] == checkFor[2][dr]) then  -- Lower/upper ghost cell in another direction too.
                     table.insert(changeIdxInDir, {d, setTo[2][d]})
                     break
                  end
               end
            end
         end
         for _, dC in ipairs(changeIdxInDir) do
            self.phiStencilIdx[nI][dC[1]] = dC[2]
         end
      end
   end

   -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil.
   for i = 1, self.phiStencilSize do
      grid:setIndex(self.phiStencilIdx[i])
      grid:getDx(self.dxBuf)
      self.dxStencil[i] = self.dxBuf:data()

      self.phiFld:fill(indexer(self.phiStencilIdx[i]), phiItr)
      self.phiStencil[i] = phiItr:data()
   end

   self._relaxation[self:idxToStencil(idx,self.cellsN)](self.dxStencil:data(), self.bcValue:data(), self.phiStencil:data(), outItr:data())

end


return ConstDiffusion
