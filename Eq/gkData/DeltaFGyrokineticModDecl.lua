-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into DeltaFGyrokinetic C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" , ["tensor"] = "Tensor"}

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, Bvars, linear)
   local emString = ""
   local geoType = "GenGeo"
   local funcSign
   funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f, double *out)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *f0, const double *f, double *out)"
   end

   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do bvarString = bvarString .. v end

   local linStr = ""
   if linear then linStr = "Linear" end

   local funcType = "double"
   local funcNm = string.format("%s%sDeltaFGyrokinetic%sVol%dx%dv%sP%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, Bvars, linear)
   local funcType  = "double"
   local geoType = "GenGeo"
   local emString  = ""
   local funcSign
   funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *fL, const double *fR, double *outL, double *outR)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR)"
   end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do bvarString = bvarString .. v end

   linStr = ""
   if linear then linStr = "Linear" end

   local funcNm = {}
   if CDIM == 1 and VDIM <= 2 then
      funcNm[1] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_x_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_vpar_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 2 then
      funcNm[1] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_x_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_y_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[3] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_vpar_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 3 and VDIM == 2 then
      funcNm[1] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_x_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_y_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[3] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_z_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[4] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_vpar_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 0 then 
      funcNm[1] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_x_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%s%sDeltaFGyrokinetic%sSurf%dx%dv%s_y_P%d", emString, linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   else
      assert(false, "DeltaFGyrokinetic equation not implemented for this dimensionality!")
   end

   local CDefStr = ""
   for d = 1, #funcNm do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, #funcNm do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end

function _M.selectStep2Vol(basisNm, CDIM, VDIM, polyOrder, geoType)
   local funcType = "double"
   linStr = ""
   if linear then linStr = "Linear" end
   local funcNm = string.format("Em%sDeltaFGyrokinetic%sStep2Vol%dx%dv%sP%d", linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f0, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectStep2Surf(basisNm, CDIM, VDIM, polyOrder, positivity, Bvars, geoType)
   local funcType  = "double"
   linStr = ""
   if linear then linStr = "Linear" end
   local funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR)"

   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do bvarString = bvarString .. v end

   local funcNm
   if polyOrder > 1 then
      funcNm = string.format("Em%sDeltaFGyrokinetic%sSurf%dx%dv%s_vpar_P%d", linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   else
      funcNm = string.format("Em%sDeltaFGyrokinetic%sSurf%dx%dv%sStep2_vpar_P%d", linStr, geoType, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathDeltaPhi(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("calcSheathDeltaPhi%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *phi, const double *phiWall, const double zVal)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathPartialReflection(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("calcSheathPartialReflectionScaled%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathReflection(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("calcSheathReflection%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
