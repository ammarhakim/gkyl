-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic Lenard-Bernstein operator equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local GkLBOModDecl = require "Eq.lboData.GkLBOModDecl"
local ffi          = require "ffi"
local xsys         = require "xsys"
local EqBase       = require "Eq.EqBase"

-- for incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Gyrokinetic Lenard-Bernstein equation on a rectangular mesh.
local GkLBO = Proto(EqBase)

-- ctor
function GkLBO:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.GkLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis  = assert(
      tbl.confBasis, "Eq.GkLBO: Must specify configuration-space basis functions to use using 'confBasis'")
   
   assert(tbl.nu, "Eq.GkLBO: Must specify collision freq using 'nu'")
   self._inNu = tbl.nu

   assert(tbl.mass, "Eq.GkLBO: Must pass mass using 'mass'")
   self._inMass = tbl.mass

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- functions to perform LBO updates.
   if self._inNu then
      self._volUpdate  = GkLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = GkLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = GkLBOModDecl.selectConstNuBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volUpdate  = GkLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = GkLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = GkLBOModDecl.selectBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- bulk velocity times collisionality field object and pointers to cell values.
   self._uNu = nil
   -- thermal speed squared times collisionality field object and pointers to cell values.
   self._vthSqNu = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._BmagInvPtr, self._BmagInvIdxr = nil, nil
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
function GkLBO:numEquations() return 1 end
function GkLBO:numWaves() return 1 end
function GkLBO:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- flux in direction dir.
function GkLBO:flux(dir, qIn, fOut)
   assert(false, "GkLBO:flux: NYI!")
end

-- Riemann problem for Gk LBO equation.
function GkLBO:rp(dir, delta, ql, qr, waves, s)
   assert(false, "GkLBO:rp: NYI!")
end

-- Compute q-fluctuations.
function GkLBO:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "GkLBO:qFluctuations: NYI!")
end

-- Maximum wave speed.
function GkLBO:maxSpeed(dir, w, dx, q)
   assert(false, "GkLBO:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function GkLBO:volTerm(w, dx, idx, q, out)
   self._BmagInv:fill(self._BmagInvIdxr(idx), self._BmagInvPtr) -- get pointer to BmagInv field.
   self._uNu:fill(self._uNuIdxr(idx), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idx), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if self._inNu then
      cflFreq = self._volUpdate(self._inMass, w:data(), dx:data(), self._BmagInvPtr:data(), self._inNu, self._uNuPtr:data(), self._vthSqNuPtr:data(), q:data(), out:data())
   else
     self._nu:fill(self._nuIdxr(idx), self._nuPtr) -- get pointer to nu field.
     cflFreq = self._volUpdate(self._inMass, w:data(), dx:data(), self._BmagInvPtr:data(), self._nuPtr:data(), self._uNuPtr:data(), self._vthSqNuPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function GkLBO:surfTerm(dir, cfl, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to uNu and vthSqNu fields.
   self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr) -- get pointer to BmagInv field.
   self._uNu:fill(self._uNuIdxr(idxl), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idxl), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if dir > self._cdim then
     if self._inNu then
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), self._BmagInvPtr:data(), self._inNu, maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     else
       self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), self._BmagInvPtr:data(), self._nuPtr:data(), maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function GkLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to uNu and vthSqNu fields.
   self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr) -- get pointer to BmagInv field.
   self._uNu:fill(self._uNuIdxr(idxl), self._uNuPtr) -- get pointer to uNu field.
   self._vthSqNu:fill(self._vthSqNuIdxr(idxl), self._vthSqNuPtr) -- get pointer to vthSqNu field.
   if dir > self._cdim then
     if self._inNu then
       vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
          self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._inNu, maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     else
       self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
       vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
          self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._nuPtr:data(), maxs, self._uNuPtr:data(), self._vthSqNuPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

function GkLBO:setAuxFields(auxFields)
   -- single aux field that has the full uNu field.
   self._BmagInv = auxFields[1]
   self._uNu     = auxFields[2]
   self._vthSqNu = auxFields[3]
   if not self._inNu then
     self._nu  = auxFields[4]
   end

   if self._isFirst then
      -- allocate pointers to field object.
      self._BmagInvPtr  = self._BmagInv:get(1)
      self._BmagInvIdxr = self._BmagInv:genIndexer()
      
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

function GkLBO:pdim() return self._pdim end
function GkLBO:cdim() return self._cdim end
function GkLBO:vdim() return self._vdim end
function GkLBO:qbym() return self._qbym end

function GkLBO:volUpdate()
   return self._volUpdate
end
function GkLBO:surfUpdate()
   return self._surfUpdate
end
function GkLBO:boundarySurfUpdate()
   return self._boundarySurfUpdate
end

return GkLBO
