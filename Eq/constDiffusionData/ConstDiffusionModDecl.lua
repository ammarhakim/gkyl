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

   local diffDirsStr = ""
   for _, d in ipairs(diffDirsIn) do diffDirsStr = diffDirsStr .. d end

   local funcType = "void"
   local funcNm   = {}
   for d = 1, #diffDirsIn do
      funcNm[d] = string.format("Const%sSurf%s%dx%sP%d_diffDirs%s_X%d", diffStr, posStr, DIM, 
                                basisNmMap[basisNm], polyOrder, diffDirsStr, diffDirsIn[d])
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

   local diffDirsStr = ""
   for _, d in ipairs(diffDirsIn) do diffDirsStr = diffDirsStr .. d end

   local funcType = "void"
   local funcNm   = {}
   for d = 1, #diffDirsIn do
      funcNm[d] = string.format("Const%sBoundarySurf%dx%sP%d_diffDirs%s_X%d", diffStr, DIM, 
                                basisNmMap[basisNm], polyOrder, diffDirsStr, diffDirsIn[d])
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

return _M
