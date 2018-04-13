-- Gkyl ------------------------------------------------------------------------
--
-- Constant Diffusion equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local ConstDiffusionModDecl = require "Eq.ConstDiffusionData.ConstDiffusionModDecl"
local ffi = require "ffi"
local xsys = require "xsys"

-- ConstDiffusion equation on a rectangular mesh
local ConstDiffusion = Proto()

-- ctor
function ConstDiffusion:init(tbl)
   
   self._nu = assert(tbl.nu, "Eq.ConstDiffusion: must specify diffusion coefficient using 'nu' ")

   -- if specified, store basis functions
   self._basis = assert(
      tbl.basis, "Eq.ConstDiffusion: Must specify basis functions to use using 'basis'")

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
   if q[0] > 0.0 then
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
   local cflFreqDiffusion = self._volUpdate(w:data(), dx:data(), self._nu, q:data(), out:data())
   return cflFreqDiffusion
end

-- Surface integral term for use in DG scheme
function ConstDiffusion:surfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._surfUpdate[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nu, ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end

function ConstDiffusion:dim() return self._ndim end

function ConstDiffusion:volUpdate()
   return self._volUpdate
end
function ConstDiffusion:surfUpdate()
   return self._surfUpdate
end

return ConstDiffusion
