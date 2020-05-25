-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into multigrid Poisson solver C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lume = require "Lib.lume"
local ffi  = require "ffi"
local _    = require "Updater.mgPoissonCalcData._MGpoissonCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

local dirLabelsLC = {'x', 'y', 'z'}
local dirLabelsUC = {'X', 'Y', 'Z'}
local boundLabel  = {'L', 'U'}

local function getStencilStrs(dimIn, bcKinds, isDG, intUpOnly) 
   -- Create a table with the strings that identify each kind of stencil location.
   -- This function assumes bcKinds is a table with dimIn entries, each one
   -- being a 2-element table where the first element is the lower boundary
   -- BC type, and the second element the upper boundary BC type.

   local stencilCount = 1
   local stencilStrs  = {}

   -- Translate the BC ID numbers to a string. For periodic (=0) we
   -- will simply use an interior stencil, assuming that the ghost cells
   -- are filled appropriately. For Dirichlet, Neumann and Robin we will use
   -- Robin kernels with appropriately chosen parameters.
   local translateBCid
   if isDG then
      translateBCid = {[0] = "", [1] = "Robin", [2] = "Robin", [3] = "Robin"}
   else
      -- FEM currently has different kernels for Dirichlet/Neumann.
      translateBCid = {[0] = "", [1] = "Dirichlet", [2] = "Neumann", [3] = "Robin", [9] = "NonPeriodic"}
   end

   local width, loOff, upOff
   if intUpOnly then
      -- Only return strings for interior and upper boundary stencils.
      width = 2
      loOff = 1
      upOff = 0
   else
      width = 3
      loOff = 0
      upOff = 0
   end

   stencilStrs[1] = ""
   for d = 1, dimIn do 
      for prevK = 1, width^(d-1) do
         for bI = 1+loOff,2-upOff do
            stencilCount = stencilCount+1
            if bcKinds[d][bI] == 0 then
               -- If BC is periodic, don't add another string. Use interior stencil.
               stencilStrs[stencilCount] = stencilStrs[prevK]
            else
               stencilStrs[stencilCount] = stencilStrs[prevK] .. boundLabel[bI] .. dirLabelsLC[d] .. translateBCid[bcKinds[d][bI]]
            end
         end
      end
   end
   for sI = 1, (width)^dimIn do   -- Append a final underscore if there are non-periodic BCs.
      if (#stencilStrs[sI] > 0) then stencilStrs[sI]=stencilStrs[sI] .. "_" end 
   end
   return stencilStrs
end

local function interiorNonPeriodicBCs(bcIn)
  -- Given the BCs in bcIn, return a new table that only
  -- differentiates between interior and non-periodic cells.
  local intUpBCs = lume.deepclone(bcIn)
  for d = 1, #bcIn do for bI = 1,2 do if intUpBCs[d][bI] ~= 0 then intUpBCs[d][bI]=9 end end end
  return intUpBCs
end

-- Select restriction operator kernel.
function _M.selectRestriction(solverKind, basisNm, dim, polyOrder, bcTypes, isDG)
   local restrictKernels = {}
   if isDG then
      local tmp = ffi.C[string.format("MGpoisson%sRestrict%dx%s_P%d", solverKind, dim, basisNmMap[basisNm], polyOrder)]
      restrictKernels[1] = tmp
   else
      local restrictStencilStr = getStencilStrs(dim, bcTypes, isDG)
      for sI = 1, 3^dim do
         local tmp = ffi.C[string.format("MGpoisson%sRestrict%dx%s_%sP%d", solverKind, dim, basisNmMap[basisNm], restrictStencilStr[sI], polyOrder)]
         restrictKernels[sI] = tmp
      end
   end
   return restrictKernels
end

-- Select prolongation operator kernel.
function _M.selectProlongation(solverKind, basisNm, dim, polyOrder, bcTypes, isDG)
   local prolongKernels = {}
   if isDG then
      local tmp = ffi.C[string.format("MGpoisson%sProlong%dx%s_P%d", solverKind, dim, basisNmMap[basisNm], polyOrder)]
      prolongKernels[1] = tmp
   else
      local prolongStencilStr = getStencilStrs(dim, bcTypes, isDG)
      for sI = 1, 3^dim do
         local tmp = ffi.C[string.format("MGpoisson%sProlong%dx%s_%sP%d", solverKind, dim, basisNmMap[basisNm], prolongStencilStr[sI], polyOrder)]
         prolongKernels[sI] = tmp
      end
   end
   return prolongKernels
end

-- Select DG to FEM kernels.
function _M.selectDGtoFEM(basisNm, dim, polyOrder, bcTypes)
   local dgToFEMkernels = {}
   local newBCs         = interiorNonPeriodicBCs(bcTypes)
   -- Create a 3^dim hypertable to place interior and upper boundary kernels.
   local dgToFEMstencilStr = getStencilStrs(dim, newBCs, false)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoissonFEM_DGtoFEM_%dx%s_%sP%d", dim, basisNmMap[basisNm], dgToFEMstencilStr[sI], polyOrder)]
      dgToFEMkernels[sI] = tmp
   end
   return dgToFEMkernels
end

-- Select restriction operator kernel.
function _M.selectFEMprojection(basisNm, dim, polyOrder, bcTypes)
   local projectKernels = {}
   local newBCs         = interiorNonPeriodicBCs(bcTypes)
   local projStencilStr = getStencilStrs(dim, newBCs, false)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoissonFEMproject%dx%s_%sP%d", dim, basisNmMap[basisNm], projStencilStr[sI], polyOrder)]
      projectKernels[sI] = tmp
   end
   return projectKernels
end

-- Select relaxation kernels.
function _M.selectRelaxation(solverKind, basisNm, dim, polyOrder, kindOfRelax, bcTypes, isDG)
   local relaxKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local relaxStencilStr = getStencilStrs(dim, bcTypes, isDG)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoisson%s%s%dx%s_%sP%d", solverKind, kindOfRelax, dim, basisNmMap[basisNm], relaxStencilStr[sI], polyOrder)]
      relaxKernels[sI] = tmp
   end
   return relaxKernels
end

-- Select residue kernels.
function _M.selectResidueCalc(solverKind, basisNm, dim, polyOrder, bcTypes, isDG)
   local residueKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local resStencilStr = getStencilStrs(dim, bcTypes, isDG)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoisson%sResidue%dx%s_%sP%d", solverKind, dim, basisNmMap[basisNm], resStencilStr[sI], polyOrder)]
      residueKernels[sI] = tmp
   end
   return residueKernels
end

-- Select kernels to compute norms.
function _M.selectFEMnorm(normType, basisNm, dim, polyOrder, bcTypes)
   local normKernels = {}
   local newBCs      = interiorNonPeriodicBCs(bcTypes)
   -- Create a 2^dim hypertable to place interior (also for lower boundary) and upper boundary kernels.
   local normStencilStr = getStencilStrs(dim, newBCs, false, true)
   for sI = 1, 2^dim do
      local tmp = ffi.C[string.format("MGpoissonFEM%snorm%dx%s_%sP%d", normType, dim, basisNmMap[basisNm], normStencilStr[sI], polyOrder)]
      normKernels[sI] = tmp
   end
   return normKernels
end

-- Select kernel that accumulates a field and a const.
function _M.selectAccuConst(basisNm, dim, polyOrder, bcTypes, isDG)
   local accuKernels = {}
   if isDG then
      local tmp = ffi.C[string.format("MGpoissonDGaccuConst%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)]
      accuKernels[1] = tmp
   else
      local newBCs = interiorNonPeriodicBCs(bcTypes)
      -- Create a 2^dim hypertable to place interior (also for lower boundary) and upper boundary kernels.
      local normStencilStr = getStencilStrs(dim, newBCs, false, true)
      for sI = 1, 2^dim do
         local tmp = ffi.C[string.format("MGpoissonFEMaccuConst%dx%s_%sP%d", dim, basisNmMap[basisNm], normStencilStr[sI], polyOrder)]
         accuKernels[sI] = tmp
      end
   end
   return accuKernels
end

return _M
