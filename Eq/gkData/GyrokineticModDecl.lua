-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Gyrokinetic C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, Bvars)
   local emString = ""
   local funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *f, double *out)"
   end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end
   local funcType = "double"
   local funcNm = string.format("%sGyrokineticVol%dx%dv%sP%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, positivity, Bvars)
   local funcType = "double"
   local emString = ""
   local funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr)"
   if isElectromagnetic then
      emString = "Em"
      funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR)"
   end
   local posString = ""
   if positivity then
      posString = "Positivity"
      if isElectromagnetic then
         funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr)"
      end
   end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end

   local funcNm = {}
   if CDIM == 1 and VDIM <= 2 then
      funcNm[1] = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 2 then
      funcNm[1] = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[3] = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 3 and VDIM == 2 then
      funcNm[1] = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[3] = string.format("%sGyrokineticSurf%s%dx%dv%s_Z_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[4] = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   elseif CDIM == 2 and VDIM == 0 then 
      funcNm[1] = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
      funcNm[2] = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
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
   local funcType = "double"
   local funcNm = string.format("EmGyrokineticStep2Vol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectStep2Surf(basisNm, CDIM, VDIM, polyOrder, positivity, Bvars)
   local funcType = "double"
   local funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR)"
   local emString = "Em"
   local posString = ""
   if positivity then
      posString = "Positivity"
      funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr)"
   end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end
   local funcNm
   if polyOrder > 1 then
      funcNm = string.format("EmGyrokineticSurf%s%dx%dv%s_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
   else
      funcNm = string.format("EmGyrokineticSurf%s%dx%dv%sStep2_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder) .. bvarString
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
