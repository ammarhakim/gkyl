-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic Lenard-Bernstein operator equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local GkLBOModDecl = require "Eq.lboData.GkLBOModDecl"
local ffi          = require "ffi"
local xsys         = require "xsys"
local EqBase       = require "Eq.EqBase"

-- For incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Gyrokinetic Lenard-Bernstein equation on a rectangular mesh.
local GkLBO = Proto(EqBase)

-- ctor.
function GkLBO:init(tbl)

   self._phaseBasis = assert(tbl.phaseBasis, 
      "Eq.GkLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis  = assert(tbl.confBasis, 
      "Eq.GkLBO: Must specify configuration-space basis functions to use using 'confBasis'")
   
   self._vParMax   = assert(tbl.vParUpper, 
      "Eq.GkLBO: Must specify maximum velocity of vPar grid in 'vParUpper'")
   self._vParMaxSq = self._vParMax^2

   assert(tbl.mass, "Eq.GkLBO: Must pass mass using 'mass'")
   self._inMass = tbl.mass

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   local nuCellConst = tbl.useCellAverageNu
   if nuCellConst == nil then
      self._varNu       = false    -- Collisionality is spatially constant.
      self._cellConstNu = true     -- Cell-wise constant nu?
   else
      self._cellConstNu = nuCellConst    -- Specify 'useCellAverageNu=true/false' for using cellwise
                                         -- constant/expanded spatially varying collisionality.
      self._varNu               = true    -- Spatially varying nu.
      self._nuPtr, self._nuIdxr = nil, nil
   end
   if self._cellConstNu then    -- Not varying within the cell.
      self._inNu = Lin.Vec(1)
   end

   -- To obtain the cell average, multiply the zeroth coefficient by this factor.
   self._cellAvFac = 1.0/math.sqrt(2.0^self._cdim)

   -- Functions to perform LBO updates.
   if self._cellConstNu then
      self._volUpdate  = GkLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = GkLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = GkLBOModDecl.selectConstNuBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volUpdate  = GkLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = GkLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = GkLBOModDecl.selectBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- Collisionality and inverse magnetic field amplitude passed as auxiliary fields.
   self._nu      = nil
   self._BmagInv = nil
   -- Bulk velocity field object and pointers to cell values.
   self._u       = nil
   -- Thermal speed squared field object and pointers to cell values.
   self._vthSq   = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._BmagInvPtr, self._BmagInvIdxr = nil, nil
   self._uPtr, self._uIdxr             = nil, nil
   self._vthSqPtr, self._vthSqIdxr     = nil, nil

   -- Flag to indicate if we are being called for first time.
   self._isFirst = true

   self.primMomCrossLimit = 0.0
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
   self._BmagInv:fill(self._BmagInvIdxr(idx), self._BmagInvPtr) -- Get pointer to BmagInv field.
   self._u:fill(self._uIdxr(idx), self._uPtr)                   -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idx), self._vthSqPtr)       -- Get pointer to vthSq field.
   if self._cellConstNu then
      if self._varNu then
         self._nu:fill(self._nuIdxr(idx), self._nuPtr)          -- Get pointer to nu field.
         self._inNu = self._nuPtr[1]*self._cellAvFac
      end
      -- If mean flow and thermal speeds are too high or if thermal
      -- speed is negative turn the LBO off (do not call kernels).
      -- Cell average values of uPar and vthSq (mind normalization).
      local uPar0   = self._uPtr[1]*self._cellAvFac
      local vthSq0  = self._vthSqPtr[1]*self._cellAvFac
      if ((math.abs(uPar0)<self._vParMax) and
          (vthSq0>0) and (vthSq0<self._vParMaxSq)) then
         cflFreq = self._volUpdate(self._inMass, w:data(), dx:data(), self._BmagInvPtr:data(), self._inNu, self._uPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
      else
         cflFreq = 0.0
         self.primMomCrossLimit = self.primMomCrossLimit+1
      end
   else
      self._nu:fill(self._nuIdxr(idx), self._nuPtr)             -- Get pointer to nu field.
      cflFreq = self._volUpdate(self._inMass, w:data(), dx:data(), self._BmagInvPtr:data(), self._nuPtr:data(), self._uPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function GkLBO:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to u and vthSq fields.
   self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr) -- Get pointer to BmagInv field.
   self._u:fill(self._uIdxr(idxl), self._uPtr)                   -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idxl), self._vthSqPtr)       -- Get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nu:fill(self._nuIdxr(idxl), self._nuPtr)       -- Get pointer to nu field.
            self._inNu = self._nuPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         -- Cell average values of uPar and vthSq (mind normalization).
         local uPar0   = self._uPtr[1]*self._cellAvFac
         local vthSq0  = self._vthSqPtr[1]*self._cellAvFac
         if ((math.abs(uPar0)<self._vParMax) and
             (vthSq0>0) and (vthSq0<self._vParMaxSq)) then
            vMuMidMax = self._surfUpdate[dir-self._cdim](
               self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), self._BmagInvPtr:data(), self._inNu, maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nu:fill(self._nuIdxr(idxl), self._nuPtr)          -- Get pointer to nu field.
         vMuMidMax = self._surfUpdate[dir-self._cdim](
            self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), self._BmagInvPtr:data(), self._nuPtr:data(), maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function GkLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- set pointer to u and vthSq fields.
   self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr) -- Get pointer to BmagInv field.
   self._u:fill(self._uIdxr(idxl), self._uPtr)             -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idxl), self._vthSqPtr) -- Get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- Get pointer to nu field.
            self._inNu = self._nuPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         -- Cell average values of uPar and vthSq (mind normalization).
         local uPar0   = self._uPtr[1]*self._cellAvFac
         local vthSq0  = self._vthSqPtr[1]*self._cellAvFac
         if ((math.abs(uPar0)<self._vParMax) and
             (vthSq0>0) and (vthSq0<self._vParMaxSq)) then
            vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
               self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._inNu, maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- get pointer to nu field.
         vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
            self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._nuPtr:data(), maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return vMuMidMax
end

function GkLBO:setAuxFields(auxFields)
   -- Single aux field that has the full BmagInv, u, vthSq (and nu) fields.
   self._BmagInv = auxFields[1]
   self._u       = auxFields[2]
   self._vthSq   = auxFields[3]
   self._nu      = auxFields[4]

   if self._isFirst then
      -- Allocate pointers to field object.
      self._BmagInvPtr  = self._BmagInv:get(1)
      self._BmagInvIdxr = self._BmagInv:genIndexer()
      
      self._uPtr  = self._u:get(1)
      self._uIdxr = self._u:genIndexer()

      self._vthSqPtr  = self._vthSq:get(1)
      self._vthSqIdxr = self._vthSq:genIndexer()

      if self._varNu then
         self._nuPtr  = self._nu:get(1)
         self._nuIdxr = self._nu:genIndexer()
      else
         self._inNu   = self._nu
      end

      self._isFirst = false -- No longer first time.
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
