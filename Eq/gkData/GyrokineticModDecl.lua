-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Gyrokinetic C++ kernel functions based on CDIM
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
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic)
   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local emString = ""
   local funcSign
   funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out)"
   end

   local bvarString = "_Bvarsxz"

   local funcType = "double"
   local funcNm = string.format("%sGyrokineticGenGeoVol%dx%dv%sP%d", emString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, positivity)
   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local funcType  = "double"
   local emString  = ""
   local posString = ""
   local funcSign
   funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR)"
   end
   if positivity then
      posString = "Positivity"
      if isElectromagnetic then
         funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR)"
      end
   end
   local bvarString = "_Bvarsxz"

   local funcNm = {}
   if CDIM == 1 and VDIM <= 2 then
      funcNm[1] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_x_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 2 then
      funcNm[1] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_x_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[3] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   elseif CDIM == 3 and VDIM == 2 then
      funcNm[1] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_x_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[3] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_z_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[4] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 0 then 
      funcNm[1] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_x_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticGenGeoSurf%s%dx%dv%s_y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   else
      assert(false, "Gyrokinetic equation not implemented for this dimensionality!")
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

function _M.selectStep2Vol(basisNm, CDIM, VDIM, polyOrder)

   local pOrder = CDIM+VDIM==5 and 1 or polyOrder

   local funcType = "double"
   local funcNm = string.format("EmGyrokineticGenGeoStep2Vol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], pOrder)
   local funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectStep2Surf(basisNm, CDIM, VDIM, polyOrder, positivity)

   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local funcType  = "double"
   local emString  = "Em"
   local posString = ""
   local funcSign
   funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR)"
   if positivity then
      posString = "Positivity"
      funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR)"
   end

   local bvarString = "_Bvarsxz"

   local funcNm
   if pOrder > 1 then
      funcNm = string.format("EmGyrokineticGenGeoSurf%s%dx%dv%s_vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   else
      funcNm = string.format("EmGyrokineticGenGeoSurf%s%dx%dv%sStep2_vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], pOrder) .. bvarString
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathDeltaPhi(basisNm, CDIM, polyOrder)
   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local funcType = "double"
   local funcNm = string.format("calcSheathDeltaPhi%dx%s_P%d", CDIM, basisNmMap[basisNm], pOrder)
   local funcSign = "(const double *phi, const double *phiWall, const double zVal)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathPartialReflection(basisNm, CDIM, VDIM, polyOrder)

   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local funcType = "void"
   local funcNm = string.format("calcSheathPartialReflectionScaled%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], pOrder)
   local funcSign = "(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSheathReflection(basisNm, CDIM, VDIM, polyOrder)
   local pOrder
   if CDIM+VDIM==5 and polyOrder>1 then
      pOrder = 1
      print("GyrokineticModDecl: Forcing selection of 3x2v p=1 gyrokinetic kernels (instead of p=2). Not used by g0.")
   else
      pOrder = polyOrder
   end

   local funcType = "void"
   local funcNm = string.format("calcSheathReflection%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], pOrder)
   local funcSign = "(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
