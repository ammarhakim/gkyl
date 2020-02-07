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

-- Translate the BC ID numbers to a string. For periodic (=0) we
-- will simply use an interior stencil, assuming that the ghost cells
-- are filled appropriately.
local translateBCid = {[0] = "", [1] = "Dirichlet", [2] = "Neumann"}

-- Select restriction operator kernel.
function _M.selectRestriction(basisNm, dim, polyOrder)
   local funcNm = string.format("MGpoissonRestrict%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select prolongation operator kernel.
function _M.selectProlongation(basisNm, dim, polyOrder)
   local funcNm = string.format("MGpoissonProlong%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

local function getStencilStrs(dimIn, bcKinds) 
   -- Create a table with the strings that identify each kind of stencil location.
   -- This function assumes bcKinds is a table with dimIn entries, each one
   -- being a 2-element table where the first element is the lower boundary
   -- BC type, and the second element the upper boundary BC type.
   local stencilCount = 1
   local stencilStrs  = {}
   stencilStrs[1]     = ""
   for d = 1, dimIn do 
      for prevK = 1, 3^(d-1) do
         stencilCount = stencilCount+1
         if bcKinds[d][1] == 0 then
            -- If BC is periodic, don't add another string. Use interior stencil.
            stencilStrs[stencilCount] = stencilStrs[prevK]
         else
            stencilStrs[stencilCount] = stencilStrs[prevK] .. "L" .. dirLabelsLC[d] .. translateBCid[bcKinds[d][1]]
         end

         stencilCount = stencilCount+1
         if bcKinds[d][2] == 0 then 
            -- If BC is periodic, don't add another string. Use interior stencil.
            stencilStrs[stencilCount] = stencilStrs[prevK]
         else
            stencilStrs[stencilCount] = stencilStrs[prevK] .. "U" .. dirLabelsLC[d] .. translateBCid[bcKinds[d][2]]
         end
      end
   end
   for sI = 1, 3^dimIn do   -- Append a final underscore if there are non-periodic BCs.
      if (#stencilStrs[sI] > 0) then stencilStrs[sI]=stencilStrs[sI] .. "_" end 
   end
   return stencilStrs
end

-- Select relaxation kernels.
function _M.selectRelaxation(basisNm, dim, polyOrder, kindOfRelax, bcTypes)
   local relaxKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local relaxStencilStr = getStencilStrs(dim, bcTypes)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoisson%s%dx%s_%sP%d", kindOfRelax, dim, basisNmMap[basisNm], relaxStencilStr[sI], polyOrder)]
      relaxKernels[sI] = tmp
   end
   return relaxKernels
end

-- Select residue kernels.
function _M.selectResidueCalc(basisNm, dim, polyOrder, bcTypes)
   local residueKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local resStencilStr = getStencilStrs(dim, bcTypes)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("MGpoissonResidue%dx%s_%sP%d", dim, basisNmMap[basisNm], resStencilStr[sI], polyOrder)]
      residueKernels[sI] = tmp
   end
   return residueKernels
end

return _M
