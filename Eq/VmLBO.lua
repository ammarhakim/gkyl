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

   self._phaseBasis    = assert(
      tbl.phaseBasis, "Eq.VmLBO: Must specify phase-space basis functions to use using 'phaseBasis'.")
   self._confBasis     = assert(
      tbl.confBasis, "Eq.VmLBO: Must specify configuration-space basis functions to use using 'confBasis'.")
   self._vMax          = assert( 
      tbl.vUpper, "Eq.VmLBO: Must specify maximum velocity of v grid in 'vUpper'.")
   local varNuIn       = tbl.varyingNu           -- Specify if collisionality varies spatially.
   local cellConstNuIn = tbl.useCellAverageNu    -- Specify whether to use cell-wise constant collisionality.
   
   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   self._cNumBasis = self._confBasis:numBasis()

   --  Store the largest velocity in the grid to compare against vtSq.
   self._vMaxSq = self._vMax[1]
   for vd = 1,self._vdim do
      if (self._vMaxSq < self._vMax[vd]) then
        self._vMaxSq = self._vMax[vd]
      end
   end
   self._vMaxSq = self._vMaxSq^2 

   -- The default is spatially constant collisionality.
   if varNuIn==nil then
      self._varNu       = false
      self._cellConstNu = true
   else
      self._varNu       = varNuIn
      if cellConstNuIn==nil then
         self._cellConstNu = true
      else
         self._cellConstNu = cellConstNuIn
      end
   end

   if self._varNu then
      self._nuSumPtr, self._nuSumIdxr = nil, nil
   end
   if self._cellConstNu then    -- Not varying within the cell.
      self._inNuSum = Lin.Vec(1)
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

   -- Collisionality passed as auxiliary field (even if just a constant number, not a field).
   self._nuSum     = nil
   -- Bulk velocity times collisionality field object and pointers to cell values.
   self._nuUSum    = nil
   -- Thermal speed squared times collisionality field object and pointers to cell values.
   self._nuVtSqSum = nil
   -- These will be set on the first call to setAuxFields() method.
   self._nuUSumPtr, self._nuUSumIdxr       = nil, nil
   self._nuVtSqSumPtr, self._nuVtSqSumIdxr = nil, nil

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

-- Code used to determine if u and vtSq crossed limits.
function VmLBO:checkPrimMomCrossings()
   local noUcrossing = true
   for vd = 1,self._vdim do
      if (math.abs(self._nuUSumPtr[(vd-1)*self._cNumBasis+1]*self._cellAvFac)>(self._vMax[vd]*self._inNuSum)) then
         noUcrossing = false
         break
      end
   end
   local nuVtSqSum0 = self._nuVtSqSumPtr[1]*self._cellAvFac
   return noUcrossing, nuVtSqSum0
end

-- Volume integral term for use in DG scheme.
function VmLBO:volTerm(w, dx, idx, q, out)
   self._nuUSum:fill(self._nuUSumIdxr(idx), self._nuUSumPtr)          -- Get pointer to u field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idx), self._nuVtSqSumPtr) -- Get pointer to vtSq field.
   if self._cellConstNu then
      if self._varNu then
         self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)    -- Get pointer to nu field.
         self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
      end
      -- If mean flow and thermal speeds are too high or if thermal
      -- speed is negative turn the LBO off (do not call kernels).
      local uCrossingNotFound, nuVtSqSum0 = self:checkPrimMomCrossings()
      if (uCrossingNotFound and (nuVtSqSum0>0) and (nuVtSqSum0<(self._vMaxSq*self._inNuSum))) then
         cflFreq = self._volUpdate(w:data(), dx:data(), self._inNuSum, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), q:data(), out:data())
      else
         cflFreq = 0.0
      end
   else
      self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)    -- Get pointer to nu field.
      cflFreq = self._volUpdate(w:data(), dx:data(), self._nuSumPtr:data(), self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function VmLBO:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- Set pointer to u and vthSq fields.
   self._nuUSum:fill(self._nuUSumIdxr(idxl), self._nuUSumPtr)          -- Get pointer to u field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idxl), self._nuVtSqSumPtr) -- Get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)    -- Get pointer to nu field.
            self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         local uCrossingNotFound, nuVtSqSum0 = self:checkPrimMomCrossings()
         if (uCrossingNotFound and (nuVtSqSum0>0) and (nuVtSqSum0<(self._vMaxSq*self._inNuSum))) then
            vMuMidMax = self._surfUpdate[dir-self._cdim](
               wl:data(), wr:data(), dxl:data(), dxr:data(), self._inNuSum, maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)    -- Get pointer to nu field.
         vMuMidMax = self._surfUpdate[dir-self._cdim](
            wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuSumPtr:data(), maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function VmLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local vMuMidMax = 0.0
   -- Set pointer to u and vthSq fields.
   self._nuUSum:fill(self._nuUSumIdxr(idxl), self._nuUSumPtr)          -- Get pointer to u field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idxl), self._nuVtSqSumPtr) -- Get pointer to vthSq field.
   if dir > self._cdim then
      if self._cellConstNu then
         if self._varNu then
            self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr) -- Get pointer to nu field.
            self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         local uCrossingNotFound, nuVtSqSum0 = self:checkPrimMomCrossings()
         if (uCrossingNotFound and (nuVtSqSum0>0) and (nuVtSqSum0<(self._vMaxSq*self._inNuSum))) then
            vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
               wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._inNuSum, maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      else
         self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr) -- Get pointer to nu field.
         vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
            wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nuSumPtr:data(), maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
     end
   end
   return vMuMidMax
end

function VmLBO:setAuxFields(auxFields)
   -- Single aux field that has the full u, vthSq (and nu) fields.
   self._nuUSum    = auxFields[1]
   self._nuVtSqSum = auxFields[2]
   self._nuSum     = auxFields[3]

   if self._isFirst then
      -- Allocate pointers to field object.
      self._nuUSumPtr  = self._nuUSum:get(1)
      self._nuUSumIdxr = self._nuUSum:genIndexer()

      self._nuVtSqSumPtr  = self._nuVtSqSum:get(1)
      self._nuVtSqSumIdxr = self._nuVtSqSum:genIndexer()

      if self._varNu then
         self._nuSumPtr  = self._nuSum:get(1)
         self._nuSumIdxr = self._nuSum:genIndexer()
      else
         self._inNuSum   = self._nuSum
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
