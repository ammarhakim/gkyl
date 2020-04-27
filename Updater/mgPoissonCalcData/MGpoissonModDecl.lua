-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into multigrid Poisson solver C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.mgPoissonCalcData._MGpoissonCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

local dirLabelsLC = {'x', 'y', 'z'}
local dirLabelsUC = {'X', 'Y', 'Z'}
local boundLabel  = {'L', 'U'}

local function getStencilStrs(dimIn, bcKinds, isDG_FEMtranslation) 
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
   if isDG_FEMtranslation then
      translateBCid = {[0] = "", [1] = "", [2] = "", [3] = ""}
   else
      translateBCid = {[0] = "", [1] = "Robin", [2] = "Robin", [3] = "Robin"}
   end

   stencilStrs[1] = ""
   for d = 1, dimIn do 
      for prevK = 1, 3^(d-1) do
         for bI = 1,2 do
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
   for sI = 1, 3^dimIn do   -- Append a final underscore if there are non-periodic BCs.
      if (#stencilStrs[sI] > 0) then stencilStrs[sI]=stencilStrs[sI] .. "_" end 
   end
   return stencilStrs
end

-- Select restriction operator kernel.
function _M.selectRestriction(solverKind, basisNm, dim, polyOrder, bcTypes, isDG)
   local restrictKernels = {}
   if isDG then
      local tmp = string.format("MGpoisson%sRestrict%dx%s_P%d", solverKind, dim, basisNmMap[basisNm], polyOrder)
      restrictKernels[1] = tmp
   else
      restrictKernels[1] = nil   -- not available yet.
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
      local prolongStencilStr = getStencilStrs(dim, bcTypes, false)
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
   -- Create a 3^dim hypertable to place interior and upper boundary kernels.
   local dgToFEMstencilStr = getStencilStrs(dim, bcTypes, true)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoissonFEM_DGtoFEM_%dx%s_%sP%d", dim, basisNmMap[basisNm], dgToFEMstencilStr[sI], polyOrder)]
      dgToFEMkernels[sI] = tmp
   end
   return dgToFEMkernels
end

-- Select relaxation kernels.
function _M.selectRelaxation(solverKind, basisNm, dim, polyOrder, kindOfRelax, bcTypes)
   local relaxKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local relaxStencilStr = getStencilStrs(dim, bcTypes, false)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoisson%s%s%dx%s_%sP%d", solverKind, kindOfRelax, dim, basisNmMap[basisNm], relaxStencilStr[sI], polyOrder)]
      relaxKernels[sI] = tmp
   end
   return relaxKernels
end

-- Select residue kernels.
function _M.selectResidueCalc(solverKind, basisNm, dim, polyOrder, bcTypes)
   local residueKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local resStencilStr = getStencilStrs(dim, bcTypes, false)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoisson%sResidue%dx%s_%sP%d", solverKind, dim, basisNmMap[basisNm], resStencilStr[sI], polyOrder)]
      residueKernels[sI] = tmp
   end
   return residueKernels
end

return _M
