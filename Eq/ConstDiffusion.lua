-- Gkyl ------------------------------------------------------------------------
--
-- Constant Diffusion equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"
local ffi = require "ffi"
local xsys = require "xsys"

-- ConstDiffusion equation on a rectangular mesh
local ConstDiffusion = Proto()

-- ctor
function ConstDiffusion:init(tbl)
   
   -- Read diffusion coefficient vector
   self._nu = Lin.Vec(3)
   assert(tbl.Dcoeff, "Eq.constDiffusion: must specify diffusion coefficient vector using 'Dcoeff' ")
   for d = 1, #tbl.Dcoeff do
      self._nu[d] = tbl.Dcoeff[d]
   end

   -- if specified, store basis functions
   self._basis = assert(
      tbl.basis, "Eq.constDiffusion: Must specify basis functions to use using 'basis'")

   -- store pointers to C kernels implementing volume and surface
   -- terms
   self._volTerm, self._surfTerms = nil, nil
   if self._basis then
      local nm, ndim, p = self._basis:id(), self._basis:ndim(), self._basis:polyOrder()
      self._volTerm = ConstDiffusionModDecl.selectVol(nm, ndim, p)
      self._surfTerms = ConstDiffusionModDecl.selectSurf(nm, ndim, p)
   end

   -- flag to indicate if we are being called for first time
   self._isFirst = true
end

-- Methods
function ConstDiffusion:numEquations() return 1 end
function ConstDiffusion:numWaves() return 1 end
function ConstDiffusion:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- flux in direction dir
function ConstDiffusion:flux(dir, qIn, fOut)
   assert(false, "ConstDiffusion:flux: NYI!")
end

-- Riemann problem for ConstDiffusion equation
function ConstDiffusion:rp(dir, delta, ql, qr, waves, s)
   assert(false, "ConstDiffusion:rp: NYI!")
end

-- Compute q-fluctuations
function ConstDiffusion:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "ConstDiffusion:qFluctuations: NYI!")
end

-- Maximum wave speed
function ConstDiffusion:maxSpeed(dir, w, dx, q)
   assert(false, "ConstDiffusion:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme
function ConstDiffusion:volTerm(w, dx, idx, q, out)
   -- streaming term
   local cflFreqDiffusion = self._volTerm(w:data(), dx:data(), self._nu:data(), q:data(), out:data())
   return cflFreqDiffusion
end

-- Surface integral term for use in DG scheme
function ConstDiffusion:surfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end

function ConstDiffusion:setAuxFields(auxFields)
end


return ConstDiffusion
