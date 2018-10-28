-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov Lenard-Bernstein operator equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local VmLBOModDecl = require "Eq.lboData.VmLBOModDecl"
local ffi          = require "ffi"
local xsys         = require "xsys"
local EqBase       = require "Eq.EqBase"

-- For incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Vlasov Lenard-Bernstein equation on a rectangular mesh.
local VmLBO = Proto(EqBase)

-- ctor.
function VmLBO:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.VmLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis  = assert(
      tbl.confBasis, "Eq.VmLBO: Must specify configuration-space basis functions to use using 'confBasis'")
   self._vMax       = assert(tbl.vUpper, "Eq.VmLBO: Must specify maximum velocity of v grid in 'vUpper'")
   
   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   self._cNumBasis = self._confBasis:numBasis()

   --  Store the largest velocity in the grid to compare against vthSq.
   self._vMaxSq = self._vMax[1]
   for vd = 1,self._vdim do
      if (self._vMaxSq < self._vMax[vd]) then
        self._vMaxSq = self._vMax[vd]
      end
   end
   self._vMaxSq = self._vMaxSq^2 

   local constNu = tbl.nu
   if constNu then
      self._varNu       = false    -- Collisionality is spatially constant.
      self._inNu        = Lin.Vec(1)
      self._inNu        = constNu
      self._cellConstNu = true     -- Collisionality is cell-wise constant.
   else
      self._cellConstNu = assert(tbl.useCellAverageNu, "Eq.VmLBO: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self._varNu               = true    -- Spatially varying nu.
      self._nu                  = nil
      self._nuPtr, self._nuIdxr = nil, nil
      if self._cellConstNu then    -- Not varying within the cell.
         self._inNu      = Lin.Vec(1)
      end
   end

   -- To obtain the cell average, multiply the zeroth coefficient by this factor.
   self._cellAvFac = 1.0/(math.sqrt(2.0^self._cdim))

   -- functions to perform LBO updates.
   if self._cellConstNu then
      self._volUpdate  = VmLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = VmLBOModDecl.selectConstNuBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volUpdate  = VmLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate = VmLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._boundarySurfUpdate = VmLBOModDecl.selectBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- bulk velocity times collisionality field object and pointers to cell values.
   self._u     = nil
   -- thermal speed squared times collisionality field object and pointers to cell values.
   self._vthSq = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._uPtr, self._uIdxr         = nil, nil
   self._vthSqPtr, self._vthSqIdxr = nil, nil

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

-- Code used to determine if u and vthSq crossed limits.
function VmLBO:checkPrimMomCrossings()
   local noUcrossing = true
   for vd = 1,self._vdim do
      if (math.abs(self._uPtr[(vd-1)*self._cNumBasis+1]*self._cellAvFac)>self._vMax[vd]) then
         noUcrossing = false
         break
      end
   end
   local vtSq0 = self._vthSqPtr[1]*self._cellAvFac
   return noUcrossing, vtSq0
end

-- Volume integral term for use in DG scheme.
function VmLBO:volTerm(w, dx, idx, q, out)
   self._u:fill(self._uIdxr(idx), self._uPtr)             -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idx), self._vthSqPtr) -- Get pointer to vthSq field.
   if self._cellConstNu then
      if self._varNu then
         self._nu:fill(self._nuIdxr(idx), self._nuPtr)    -- Get pointer to nu field.
         self._inNu = self._nuPtr[1]*self._cellAvFac
      end
      -- If mean flow and thermal speeds are too high or if thermal
      -- speed is negative turn the LBO off (do not call kernels).
      local uCrossingNotFound, vthSq0 = self:checkPrimMomCrossings()
      if (uCrossingNotFound and (vthSq0>0) and (vthSq0<self._vMaxSq)) then
         cflFreq = self._volUpdate(w:data(), dx:data(), self._inNu, self._uPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
      else
         cflFreq = 0.0
      end
   else
      self._nu:fill(self._nuIdxr(idx), self._nuPtr)    -- Get pointer to nu field.
      cflFreq = self._volUpdate(w:data(), dx:data(), self._nuPtr:data(), self._uPtr:data(), self._vthSqPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function VmLBO:surfTerm(dir, dt, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- Set pointer to u and vthSq fields.
   self._u:fill(self._uIdxr(idxl), self._uPtr) -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idxl), self._vthSqPtr) -- Get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nu:fill(self._nuIdxr(idxl), self._nuPtr)    -- Get pointer to nu field.
            self._inNu = self._nuPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         local uCrossingNotFound, vthSq0 = self:checkPrimMomCrossings()
         if (uCrossingNotFound and (vthSq0>0) and (vthSq0<self._vMaxSq)) then
            vMuMidMax = self._surfUpdate[dir-self._cdim](
               wl:data(), wr:data(), dxl:data(), dxr:data(), self._inNu, maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nu:fill(self._nuIdxr(idxl), self._nuPtr)    -- Get pointer to nu field.
         vMuMidMax = self._surfUpdate[dir-self._cdim](
            wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuPtr:data(), maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function VmLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- Set pointer to u and vthSq fields.
   self._u:fill(self._uIdxr(idxl), self._uPtr) -- Get pointer to u field.
   self._vthSq:fill(self._vthSqIdxr(idxl), self._vthSqPtr) -- get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- Get pointer to nu field.
            self._inNu = self._nuPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         local uCrossingNotFound, vthSq0 = self:checkPrimMomCrossings()
         if (uCrossingNotFound and (vthSq0>0) and (vthSq0<self._vMaxSq)) then
            vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
               wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._inNu, maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nu:fill(self._nuIdxr(idxl), self._nuPtr) -- Get pointer to nu field.
         vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
            wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nuPtr:data(), maxs, self._uPtr:data(), self._vthSqPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

function VmLBO:setAuxFields(auxFields)
   -- Single aux field that has the full u, vthSq (and nu) fields.
   self._u     = auxFields[1]
   self._vthSq = auxFields[2]
   if self._varNu then
      self._nu = auxFields[3]
   end

   if self._isFirst then
      -- Allocate pointers to field object.
      self._uPtr  = self._u:get(1)
      self._uIdxr = self._u:genIndexer()

      self._vthSqPtr  = self._vthSq:get(1)
      self._vthSqIdxr = self._vthSq:genIndexer()

      if self._varNu then
         self._nuPtr  = self._nu:get(1)
         self._nuIdxr = self._nu:genIndexer()
      end

      self._isFirst = false -- No longer first time.
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
