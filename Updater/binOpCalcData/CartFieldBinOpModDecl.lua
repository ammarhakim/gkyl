-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into binary operator calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

ffi.cdef [[
typedef struct binOpData_t binOpData_t;
binOpData_t* new_binOpData_t(int nbasis_S, int nbasis_D);
]]

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select kernel function to compute specified operation
-- between inputs with same dimensionality.
function _M.selectBinOpCalcS(op, basisNm, CDIM, VDIM, polyOrder, applyPos)
   local funcType = "void"
   local posString = ""
   if (applyPos and op=="Divide") then posString = "Positivity" end
   local funcNm
   if VDIM then 
      funcNm = string.format("CartFieldBinOp%s%s%dx%dv%s_P%d", op, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   else
      funcNm = string.format("CartFieldBinOp%s%s%dx%s_P%d", op, posString, CDIM, basisNmMap[basisNm], polyOrder)
   end
   local funcSign = "(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select kernel function to compute specified operation
-- between inputs of different dimensionality.
function _M.selectBinOpCalcD(op, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("CartFieldBinOpConfPhase%s%dx%dv%s_P%d", op, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
