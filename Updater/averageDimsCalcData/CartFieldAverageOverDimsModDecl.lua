--
-- Dispatch into ConstDiffusion C++ kernel functions based on DIM,
-- basis functions, polyOrder and the list of diffusive directions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["tensor"] = "tensor" }

local _M = {}

-- Select kernel to compute average.
function _M.selectAverageOverDirsKernel(DIM, basisNm, polyOrder, avgDirsIn)
   local avgDirsStr = ""
   for _, d in ipairs(avgDirsIn) do avgDirsStr = avgDirsStr .. d end

   local funcType = "void"
   local funcNm   = string.format("avg_over_dims_%dx_p%d_%s_avgDirs%s", DIM, polyOrder, basisNmMap[basisNm], avgDirsStr)
   local funcSign = "(double rNumCells, const double *fIn, double *fAvgOut)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
