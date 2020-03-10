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
local DataStruct   = require "DataStruct"
local Updater      = require "Updater"

-- For incrementing in updater.
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- Gyrokinetic Lenard-Bernstein equation on a rectangular mesh.
local GkLBO = Proto(EqBase)

-- ctor.
function GkLBO:init(tbl)

   self._grid          = assert(tbl.onGrid, 
      "Eq.GkLBO: Must specify the grid using 'onGrid'.")
   self._phaseBasis    = assert(tbl.phaseBasis, 
      "Eq.GkLBO: Must specify phase-space basis functions to use using 'phaseBasis'.")
   self._confBasis     = assert(tbl.confBasis, 
      "Eq.GkLBO: Must specify configuration-space basis functions to use using 'confBasis'.")
   
   self._vParMax       = assert(tbl.vParUpper, 
      "Eq.GkLBO: Must specify maximum velocity of vPar grid in 'vParUpper'.")
   self._vParMaxSq     = self._vParMax^2
   local varNuIn       = tbl.varyingNu    -- Specify if collisionality varies spatially.
   local cellConstNuIn = tbl.useCellAverageNu    -- Specify whether to use cell-wise constant collisionality.

   assert(tbl.mass, "Eq.GkLBO: Must pass mass using 'mass'.")
   self._inMass = tbl.mass

   self._positivity = xsys.pickBool(tbl.positivity,false)   -- Positivity preserving option.

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

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
   self._cellAvFac = 1.0/math.sqrt(2.0^self._cdim)

   -- Functions to perform LBO updates.
   if self._cellConstNu then
      self._volUpdate          = GkLBOModDecl.selectConstNuVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate         = GkLBOModDecl.selectConstNuSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._positivity)
      self._boundarySurfUpdate = GkLBOModDecl.selectConstNuBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._positivity)
   else
      self._volUpdate          = GkLBOModDecl.selectVol(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfUpdate         = GkLBOModDecl.selectSurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._positivity)
      self._boundarySurfUpdate = GkLBOModDecl.selectBoundarySurf(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._positivity)
   end

   -- Collisionality and inverse magnetic field amplitude passed as auxiliary fields.
   self._nuSum     = nil
   self._BmagInv   = nil
   -- Bulk velocity field object and pointers to cell values.
   self._nuUSum    = nil
   -- Thermal speed squared field object and pointers to cell values.
   self._nuVtSqSum = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._BmagInvPtr, self._BmagInvIdxr     = nil, nil
   self._nuUSumPtr, self._nuUSumIdxr       = nil, nil
   self._nuVtSqSumPtr, self._nuVtSqSumIdxr = nil, nil

   -- Flag to indicate if we are being called for first time.
   self._isFirst = true

   self.primMomCrossLimit = 0.0

   if self._positivity then
      self.fRhsVol = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._phaseBasis:numBasis(),
         ghost         = {1, 1}
      }
      self.fRhsVol_ptr = self.fRhsVol:get(1)

      self.fRhsSurf = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._phaseBasis:numBasis(),
         ghost         = {1, 1}
      }
      self.fRhsSurf_L_ptr = self.fRhsSurf:get(1)
      self.fRhsSurf_R_ptr = self.fRhsSurf:get(1)
      self.fRhsIdxr = self.fRhsVol:genIndexer()

      self.posRescaler = Updater.PositivityRescale {
         onGrid = self._grid,
         basis  = self._phaseBasis,
      }

      self.volTermScaleFac = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = 1,   -- Only need one number per cell.
         ghost         = {1, 1}
      }
      self.volTermScaleFac_ptr = self.volTermScaleFac:get(1)
      self.volTermScaleFacIdxr = self.volTermScaleFac:genIndexer()
   end

   self.positivityWeightByDir = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._vdim + 1,
      ghost         = {1, 1},
   }
   self.positivityWeightByDirIdxr  = self.positivityWeightByDir:genIndexer()
   self.positivityWeightByDir_ptr  = self.positivityWeightByDir:get(1)
   self.positivityWeightByDirL_ptr = self.positivityWeightByDir:get(1)
   self.positivityWeightByDirR_ptr = self.positivityWeightByDir:get(1)
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
function GkLBO:volTerm(w, dx, idx, f_ptr, out_ptr)
   self._BmagInv:fill(self._BmagInvIdxr(idx), self._BmagInvPtr)          -- Get pointer to BmagInv field.
   self._nuUSum:fill(self._nuUSumIdxr(idx), self._nuUSumPtr)             -- Get pointer to sum(nu*u) field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idx), self._nuVtSqSumPtr)    -- Get pointer to sum(nu*vtSq) field.

   self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)
   if self._positivity then
      -- Positivity algorithms write to a volume term cartField stored in GkLBO.
      self.fRhsVol:fill(self.fRhsIdxr(idx), self.fRhsVol_ptr)
   else
      self.fRhsVol_ptr = out_ptr
   end

   if self._cellConstNu then
      if self._varNu then
         self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)          -- Get pointer to sum(nu) field.
         self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
      end
      -- If mean flow and thermal speeds are too high or if thermal
      -- speed is negative turn the LBO off (do not call kernels).
      -- Cell average values of uPar and vtSq (mind normalization).
      local nuUParSum0 = self._nuUSumPtr[1]*self._cellAvFac
      local nuVtSqSum0 = self._nuVtSqSumPtr[1]*self._cellAvFac
      if ((math.abs(nuUParSum0)<(self._vParMax*self._inNuSum)) and
          (nuVtSqSum0>0) and (nuVtSqSum0<(self._vParMaxSq*self._inNuSum))) then
         cflRate = self._volUpdate(self._inMass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self._BmagInvPtr:data(), self._inNuSum, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
         -- if using positivity, the volume term will be repeated in step 2, and we will use the cflFreq from there
         -- to avoid double counting, set cflFreq = 0 here
         if self._positivity then cflRate = 0.0 end
      else
         cflRate = 0.0
         self.primMomCrossLimit = self.primMomCrossLimit+1
      end
   else
      self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)             -- Get pointer to sum(nu) field.
      cflRate = self._volUpdate(self._inMass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self._BmagInvPtr:data(), self._nuSumPtr:data(), self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
   end
   return cflRate
end

-- Surface integral term for use in DG scheme.
function GkLBO:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, fL_ptr, fR_ptr, outL_ptr, outR_ptr)
   local vMuMidMax = 0.0
   if dir > self._cdim then
      -- Set pointer to BmagInv, sum(nu*u) and sum(nu*vtSq) fields.
      self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr)          -- Get pointer to BmagInv field.
      self._nuUSum:fill(self._nuUSumIdxr(idxl), self._nuUSumPtr)             -- Get pointer to sum(u) field.
      self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idxl), self._nuVtSqSumPtr)    -- Get pointer to sum(nu*vtSq) field.

      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxl), self.positivityWeightByDirL_ptr)
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxr), self.positivityWeightByDirR_ptr)

      if self._positivity then
         -- Positivity algorithms write to a surface term cartField stored in GkLBO.
         self.fRhsSurf:fill(self.fRhsIdxr(idxl), self.fRhsSurf_L_ptr)
         self.fRhsSurf:fill(self.fRhsIdxr(idxr), self.fRhsSurf_R_ptr)
      else
         self.fRhsSurf_L_ptr = outL_ptr
         self.fRhsSurf_R_ptr = outR_ptr
      end

      if self._cellConstNu then
         if self._varNu then
            self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)       -- Get pointer to sum(nu) field.
            self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
         end
         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         -- Cell average values of uPar and vtSq (mind normalization).
         local nuUParSum0 = self._nuUSumPtr[1]*self._cellAvFac
         local nuVtSqSum0 = self._nuVtSqSumPtr[1]*self._cellAvFac
         if ((math.abs(nuUParSum0)<(self._vParMax*self._inNuSum)) and
             (nuVtSqSum0>0) and (nuVtSqSum0<(self._vParMaxSq*self._inNuSum))) then
            vMuMidMax = self._surfUpdate[dir-self._cdim](
               self._inMass, self.positivityWeightByDirL_ptr:data(), self.positivityWeightByDirR_ptr:data(), wl:data(), wr:data(), dxl:data(), dxr:data(), dtApprox, self._BmagInvPtr:data(), self._inNuSum, dtApprox, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), fL_ptr:data(), fR_ptr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
         end
      else
         self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)          -- Get pointer to sum(nu) field.
         vMuMidMax = self._surfUpdate[dir-self._cdim](
            self._inMass, self.positivityWeightByDirL_ptr:data(), self.positivityWeightByDirR_ptr:data(), wl:data(), wr:data(), dxl:data(), dxr:data(), self._BmagInvPtr:data(), self._nuSumPtr:data(), dtApprox, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), fL_ptr:data(), fR_ptr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
      end
   end
   return vMuMidMax
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function GkLBO:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outL_ptr, outR_ptr)
   local vMuMidMax = 0.0
   -- Set pointer to BmagInv, sum(nu*u) and sum(nu*vtSq) fields.
   self._BmagInv:fill(self._BmagInvIdxr(idxl), self._BmagInvPtr)          -- Get pointer to BmagInv field.
   self._nuUSum:fill(self._nuUSumIdxr(idxl), self._nuUSumPtr)             -- Get pointer to sum(nu*u) field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idxl), self._nuVtSqSumPtr)    -- Get pointer to sum(nu*vtSq) field.
   if dir > self._cdim then
      if self._positivity then
         -- Positivity algorithms write to a surface term cartField stored in GkLBO.
         self.fRhsSurf:fill(self.fRhsIdxr(idxl), self.fRhsSurf_L_ptr)
         self.fRhsSurf:fill(self.fRhsIdxr(idxr), self.fRhsSurf_R_ptr)
      else
         self.fRhsSurf_L_ptr = outL_ptr
         self.fRhsSurf_R_ptr = outR_ptr
      end

      if self._cellConstNu then
         if self._varNu then
            self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)    -- Get pointer to sum(nu) field.
            self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
         end

         -- If mean flow and thermal speeds are too high or if thermal
         -- speed is negative turn the LBO off (do not call kernels).
         -- Cell average values of uPar and vtSq (mind normalization).
         local nuUParSum0 = self._nuUSumPtr[1]*self._cellAvFac
         local nuVtSqSum0 = self._nuVtSqSumPtr[1]*self._cellAvFac
         if ((math.abs(nuUParSum0)<(self._vParMax*self._inNuSum)) and
             (nuVtSqSum0>0) and (nuVtSqSum0<(self._vParMaxSq*self._inNuSum))) then
            vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
               self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._inNuSum, maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
         end
      else
         self._nuSum:fill(self._nuSumIdxr(idxl), self._nuSumPtr)    -- Get pointer to sum(nu) field.
         vMuMidMax = self._boundarySurfUpdate[dir-self._cdim](
            self._inMass, wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._BmagInvPtr:data(), self._nuSumPtr:data(), maxs, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), ql:data(), qr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
      end
   end
   return vMuMidMax
end

function GkLBO:clearRhsTerms()
   if self._positivity then
      self.fRhsVol:clear(0.0)
      self.fRhsSurf:clear(0.0)
      self.volTermScaleFac:clear(0.0)
   end
end

-- When using positivity-preserving algorithm we:
--   1) Compute the flux-limited surface terms (in surfTerm), and compute the volume terms (in volTerm).
--   2) Compute the scaling factor of the volume term (in getPositivity), but don't necessarily scale the volume term.
--   3) Compute self-prim moments again (in GkSpecies advanceStep2) with slightly different weak problem.
--   4) Compute the volume term again and apply it with new uPar and vtSq (in volTermStep2).
--   5) Rescale volume term with scaling factor from (2).
function GkLBO:getPositivityRhs(tCurr, dtApprox, fIn, fRhs)
   weightDirs = {}
   for d = 1, self._vdim do
      weightDirs[d] = d
   end
   -- fIn + fac*dt*fVol + dt*fSurf > 0.
   -- Compute the scaling factor for volume term, fac = self.volTermScaleFac, but don't do the scaling yet.
   if dtApprox > 0 and tCurr > 0 then self.posRescaler:calcVolTermRescale(tCurr, dtApprox, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurf, self.fRhsVol, self.volTermScaleFac) end
end

-- Step 2 volume integral term for use in positivity-LBO scheme.
function GkLBO:volTermStep2(w, dx, idx, f_ptr, out_ptr)
   self._BmagInv:fill(self._BmagInvIdxr(idx), self._BmagInvPtr)          -- Get pointer to BmagInv field.
   self._nuUSum:fill(self._nuUSumIdxr(idx), self._nuUSumPtr)             -- Get pointer to sum(nu*u) field.
   self._nuVtSqSum:fill(self._nuVtSqSumIdxr(idx), self._nuVtSqSumPtr)    -- Get pointer to sum(nu*vtSq) field.

   self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)

   -- Positivity algorithms write to a volume term cartField stored in GkLBO.
   self.fRhsVol:fill(self.fRhsIdxr(idx), self.fRhsVol_ptr)

   if self._cellConstNu then
      if self._varNu then
         self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)          -- Get pointer to sum(nu) field.
         self._inNuSum = self._nuSumPtr[1]*self._cellAvFac
      end
      -- If mean flow and thermal speeds are too high or if thermal
      -- speed is negative turn the LBO off (do not call kernels).
      -- Cell average values of uPar and vtSq (mind normalization).
      local nuUParSum0 = self._nuUSumPtr[1]*self._cellAvFac
      local nuVtSqSum0 = self._nuVtSqSumPtr[1]*self._cellAvFac
      if ((math.abs(nuUParSum0)<(self._vParMax*self._inNuSum)) and
          (nuVtSqSum0>0) and (nuVtSqSum0<(self._vParMaxSq*self._inNuSum))) then
         cflRate = self._volUpdate(self._inMass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self._BmagInvPtr:data(), self._inNuSum, self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
      else
         cflRate = 0.0
         self.primMomCrossLimit = self.primMomCrossLimit+1
      end
   else
      self._nuSum:fill(self._nuSumIdxr(idx), self._nuSumPtr)             -- Get pointer to sum(nu) field.
      cflRate = self._volUpdate(self._inMass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self._BmagInvPtr:data(), self._nuSumPtr:data(), self._nuUSumPtr:data(), self._nuVtSqSumPtr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
   end

   return cflRate
end

-- Step 2 surface integral term for use in positivity-LBO scheme.
function GkLBO:surfTermStep2(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, fL_ptr, fR_ptr, outL_ptr, outR_ptr)
   local vMuMidMax = 0.0
   return vMuMidMax
end

function GkLBO:getPositivityRhsStep2(tCurr, dtApprox, fIn, fRhs)
   self.fRhsVol:scaleByCell(self.volTermScaleFac)

   fRhs:accumulate(1.0, self.fRhsSurf, 1.0, self.fRhsVol)
end

function GkLBO:setPositivityWeights(cflRateByCell)
   -- set total weight = positivityWeightByDir[0] = cflRateByCell in each cell
   -- other elements will be overwritten in kernels
   self.positivityWeightByDir:clear(1.0)
   self.positivityWeightByDir:scaleByCell(cflRateByCell)
end

function GkLBO:setAuxFields(auxFields)
   -- Single aux field that has the full BmagInv, sum(nu*u), sum(nu*vtSq) and sum(nu) fields.
   self._BmagInv   = auxFields[1]
   self._nuUSum    = auxFields[2]
   self._nuVtSqSum = auxFields[3]
   self._nuSum     = auxFields[4]

   if self._isFirst then
      -- Allocate pointers to field object.
      self._BmagInvPtr  = self._BmagInv:get(1)
      self._BmagInvIdxr = self._BmagInv:genIndexer()
      
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
