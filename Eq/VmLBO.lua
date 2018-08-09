-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov Lenard-Bernstein operator equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local VmLBOModDecl = require "Eq.lboData.VmLBOModDecl"
local ffi          = require "ffi"
local xsys         = require "xsys"
local EqBase       = require "Eq.EqBase"

-- for incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Vlasov Lenard-Bernstein equation on a rectangular mesh.
local VmLBO = Proto(EqBase)

-- ctor
function VmLBO:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.VmLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis  = assert(
      tbl.confBasis, "Eq.VmLBO: Must specify configuration-space basis functions to use using 'confBasis'")
   
   assert(tbl.nu, "Eq.VmLBO: Must specify collisionality using 'nu'")
   self._inNu = Lin.Vec(1)
   self._inNu = tbl.nu

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- functions to perform LBO updates.
   if self._inNu then
      self._volUpdate  = VmLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = VmLBOModDecl.selectConstNuBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volUpdate  = VmLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = VmLBOModDecl.selectBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- bulk velocity times collisionality field object and pointers to cell values.
   self._uNu = nil
   -- thermal speed squared times collisionality field object and pointers to cell values.
   self._vthSqNu = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._uNuPtr, self._uNuIdxr = nil, nil
   self._vthSqNuPtr, self._vthSqNuIdxr = nil, nil

   if not self._inNu then
     self._nu = nil
     self._nuPtr, self._nuIdxr = nil, nil
   end

   -- flag to indicate if we are being called for first time.
   self._isFirst = true
end

-- Methods.
function VmLBO:numEquations() return 1 end
function VmLBO:numWaves() return 1 end
function VmLBO:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- flux in direction dir.
function VmLBO:flux(dir, qIn, fOut)
   assert(false, "VmLBO:flux: NYI!")
end

-- Riemann problem for Vlasov LBO equation.
function VmLBO:rp(dir, delta, ql, qr, waves, s)
   assert(false, "VmLBO:rp: NYI!")
end

-- Compute q-fluctuations.
function VmLBO:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "VmLBO:qFluctuations: NYI!")
end

-- Maximum wave speed.
function VmLBO:maxSpeed(dir, w, dx, q)
   assert(false, "VmLBO:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function VmLBO:volTerm(w, dx, idx, q, out)
   self._uNu:fill(self._uNuIdxr(idx), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idx), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if self._inNu then
     cflFreq = self._volUpdate(w:data(), dx:data(), self._inNu, self._uNuPtr:data(), self._vthSqNuPtr:data(), q:data(), out:data())
   else
     self._nu:fill(self._nuIdxr(idx), self._nuPtr) -- get pointer to nu field.
     cflFreq = self._volUpdate(w:data(), dx:data(), self._nuPtr:data(), self._uNuPtr:data(), self._vthSqNuPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function VmLBO:surfTerm(dir, dt, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to uNu and vthSqNu fields.
   self._uNu:fill(self._uNuIdxr(idxl), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idxl), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if dir > self._cdim then
     if self._inNu then
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), self._inNu, maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     else
       self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuPtr:data(), maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function VmLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to uNu and vthSqNu fields.
   self._uNu:fill(self._uNuIdxr(idxl), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idxl), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if dir > self._cdim then
     if self._inNu then
       vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._inNu, maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     else
       self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
       vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nuPtr:data(), maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

function VmLBO:setAuxFields(auxFields)
   -- single aux field that has the full uNu field.
   self._uNu     = auxFields[1]
   self._vthSqNu = auxFields[2]
   if not self._inNu then
     self._nu  = auxFields[3]
   end

   if self._isFirst then
      -- allocate pointers to field object.
      self._uNuPtr  = self._uNu:get(1)
      self._uNuIdxr = self._uNu:genIndexer()

      self._vthSqNuPtr  = self._vthSqNu:get(1)
      self._vthSqNuIdxr = self._vthSqNu:genIndexer()

      if not self._inNu then
         self._nuPtr     = self._nu:get(1)
         self._nuIdxr    = self._nu:genIndexer()
      end

      self._isFirst    = false -- no longer first time.
   end
end

function VmLBO:pdim() return self._pdim end
function VmLBO:cdim() return self._cdim end
function VmLBO:vdim() return self._vdim end
function VmLBO:qbym() return self._qbym end

function VmLBO:volUpdate()
   return self._volUpdate
end
function VmLBO:surfUpdate()
   return self._surfUpdate
end
function VmLBO:boundarySurfUpdate()
   return self._boundarySurfUpdate
end

return VmLBO
