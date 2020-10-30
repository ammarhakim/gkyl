-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into ConstDiffusion C++ kernel functions based on DIM,
-- basis functions, polyOrder and the list of diffusive directions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, DIM, polyOrder, diffDirsIn, diffOrder)
   diffStr = "Diffusion"
   if diffOrder > 2 then diffStr = "HyperDiffusion" .. diffOrder end

   local diffDirsStr = ""
   for _, d in ipairs(diffDirsIn) do diffDirsStr = diffDirsStr .. d end

   local funcType = "double"
   local funcNm   = string.format("Const%sVol%dx%sP%d_diffDirs%s", diffStr, DIM, basisNmMap[basisNm], polyOrder, diffDirsStr)
   local funcSign = "(const double *w, const double *dx, const double *nu, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, DIM, polyOrder, diffDirsIn, diffOrder, applyPos)
   diffStr = "Diffusion"
   if diffOrder > 2 then diffStr = "HyperDiffusion" .. diffOrder end

   local posStr = ""
   if applyPos then posStr = "Positivity" end

   local funcType = "void"
   local funcNm   = {}
   for d = 1, #diffDirsIn do
      funcNm[d] = string.format("Const%sSurf%s%dx%sP%d_X%d", diffStr, posStr, DIM, 
                                basisNmMap[basisNm], polyOrder, diffDirsIn[d])
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, #diffDirsIn do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, DIM do kernels[d] = nil end
   for d = 1, #diffDirsIn do
      local tmp = ffi.C[funcNm[d]]
      kernels[diffDirsIn[d]] = tmp
   end
   return kernels
end

-- Select functions to compute boundary surface terms (output is a table of functions).
function _M.selectBoundarySurf(basisNm, DIM, polyOrder, diffDirsIn, diffOrder, applyPos)
   diffStr = "Diffusion"
   if diffOrder > 2 then diffStr = "HyperDiffusion" .. diffOrder end

   local funcType = "void"
   local funcNm   = {}
   for d = 1, #diffDirsIn do
      funcNm[d] = string.format("Const%sBoundarySurf%dx%sP%d_X%d", diffStr, DIM, 
                                basisNmMap[basisNm], polyOrder, diffDirsIn[d])
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, #diffDirsIn do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, DIM do kernels[d] = nil end
   for d = 1, #diffDirsIn do
      local tmp = ffi.C[funcNm[d]]
      kernels[diffDirsIn[d]] = tmp
   end
   return kernels
end

-- Select kernels that implement BCs specific to constDiffusion term (output is a table of functions).
function _M.selectBCs(basisNm, DIM, polyOrder, diffOrder, bcType)
   diffStr = "Diffusion"
   if diffOrder > 2 then diffStr = "HyperDiffusion" .. diffOrder end

   local funcType = "void"
   local funcNm   = {}
   for d = 1, DIM do
      funcNm[d] = {string.format("Const%sBC%dx%sP%d_%s_X%dlower", diffStr, DIM, 
                                 basisNmMap[basisNm], polyOrder, bcType, d),
                   string.format("Const%sBC%dx%sP%d_%s_X%dupper", diffStr, DIM, 
                                 basisNmMap[basisNm], polyOrder, bcType, d)}
   end
   local funcSign = "(const double dx, const double *fSkin, const double fBC, double *fGhost)"

   local CDefStr = ""
   for d = 1, DIM do
      for bI = 1,2 do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d][bI] .. funcSign .. ";\n") end
   end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, DIM do kernels[d] = nil end
   for d = 1, DIM do
      local tmp  = {ffi.C[funcNm[d][1]], ffi.C[funcNm[d][2]]}
      kernels[d] = tmp
   end
   return kernels
end

local function getStencilStrs(dimIn, bcKinds)
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
   -- FEM currently has different kernels for Dirichlet/Neumann.
   translateBCid = {[0] = "", [1] = "Dirichlet", [2] = "Neumann", [3] = "Robin", [9] = "NonPeriodic"}

   local width, loOff, upOff
   width, loOff, upOff = 3, 0, 0

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

-- Select FEM diffusion kernels.
function _M.selectFEMdiff(basisNm, dim, polyOrder, bcTypes)
   local diffKernels = {}
   -- Create a 3^dim hypertable to place lower boundary, interior and upper boundary kernels.
   local diffStencilStr = getStencilStrs(dim, bcTypes)
   for sI = 1, 3^dim do
      local tmp = ffi.C[string.format("ConstDiffusionFEM%dx%s_%sP%d", dim, basisNmMap[basisNm], diffStencilStr[sI], polyOrder)]
      diffKernels[sI] = tmp
   end
   return diffKernels
end

return _M
