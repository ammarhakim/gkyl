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

-- for incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Vlasov Lenard-Bernstein equation on a rectangular mesh.
local VmLBO = Proto()

-- ctor
function VmLBO:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.VmLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Eq.VmLBO: Must specify configuration-space basis functions to use using 'confBasis'")
   
   assert(tbl.nu, "Eq.VmLBO: Must specify configuration-space basis functions to use using 'nu'")
   self._inNu = Lin.Vec(1)
   self._inNu = tbl.nu

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- functions to perform LBO updates.
   if self._inNu then
      self._volUpdate = VmLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volUpdate = VmLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- flow speed field object and pointers to cell values.
   self._flowU = nil
   -- thermal speed squared field object and pointers to cell values.
   self._vthSq = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._flowUPtr, self._flowUIdxr = nil, nil
   self._vthSqPtr, self._vthSqIdxr = nil, nil

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
   self._flowU:fill(self._flowUIdxr(idx), self._flowUPtr) -- get pointer to flowU field.
   self._vthSq:fill(self._vthSqIdxr(idx), self._vthSqPtr) -- get pointer to vthSq field.
   if self._inNu then
     cflFreq = self._volUpdate(w:data(), dx:data(), self._inNu, self._flowUPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
   else
     self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
     cflFreq = self._volUpdate(w:data(), dx:data(), self._nuPtr:data(), self._flowUPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function VmLBO:surfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to flowU and vthSq fields.
   self._flowU:fill(self._flowUIdxr(idxl), self._flowUPtr) -- get pointer to flowU field.
   self._vthSq:fill(self._vthSqIdxr(idxl), self._vthSqPtr) -- get pointer to vthSq field.
   if dir > self._cdim then
     if self._inNu then
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), self._inNu, maxs, self._flowUPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     else
       self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
       vMuMidMax = self._surfUpdate[dir-self._cdim](
          wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuPtr:data(), maxs, self._flowUPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

function VmLBO:setAuxFields(auxFields)
   -- single aux field that has the full flowU field.
   self._flowU = auxFields[1]
   self._vthSq = auxFields[2]
   if not self._inNu then
     self._nu  = auxFields[3]
   end

   if self._isFirst then
      -- allocate pointers to field object.
      self._flowUPtr  = self._flowU:get(1)
      self._flowUIdxr = self._flowU:genIndexer()

      self._vthSqPtr  = self._vthSq:get(1)
      self._vthSqIdxr = self._vthSq:genIndexer()

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

return VmLBO
