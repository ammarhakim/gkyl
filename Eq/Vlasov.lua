-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local VlasovModDecl = require "Eq.vlasovData.VlasovModDecl"
local ffi = require "ffi"
local xsys = require "xsys"

-- for incrementing in updater
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- compute emOut = q/m*emIn
local function rescaleEmField(qbym, emIn, emOut)
   ffi.C.vlasovIncr(#emOut, emIn:data(), qbym, emOut:data())
end

-- Vlasov equation on a rectangular mesh
local Vlasov = Proto()

-- ctor
function Vlasov:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.Vlasov: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Eq.Vlasov: Must specify configuration-space basis functions to use using 'confBasis'")
   
   local charge = assert(tbl.charge, "Eq.Vlasov: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Eq.Vlasov: must specify mass using 'mass' ")

   self._qbym = charge/mass -- only q/m ratio is ever needed

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- functions to perform streaming updates
   self._volStreamUpdate = VlasovModDecl.selectVolStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovModDecl.selectSurfStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   -- check if we have a electric and magnetic field
   local hasElcField = xsys.pickBool(tbl.hasElectricField, true)
   local hasMagField = xsys.pickBool(tbl.hasMagneticField, true)

   self._hasForceTerm = false -- flag to indicate if we have any force terms at all
   if hasElcField or hasMagField then
      self._hasForceTerm = true
   end

   self._volForceUpdate, self._surfForceUpdate = nil, nil
   -- functions to perform force updates from electric field
   if hasMagField then 
      self._volForceUpdate = VlasovModDecl.selectVolElcMag(
	 self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfForceUpdate = VlasovModDecl.selectSurfElcMag(
	 self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volForceUpdate = VlasovModDecl.selectVolElc(
	 self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfForceUpdate = VlasovModDecl.selectSurfElc(
	 self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- EM field object and pointers to cell values
   self._emField = nil
   -- (these will be set on the first call to setAuxFields() method)
   self._emPtr, self._emIdxr = nil, nil
   self._emAccel = nil

   -- flag to indicate if we are being called for first time
   self._isFirst = true
end

-- Methods
function Vlasov:numEquations() return 1 end
function Vlasov:numWaves() return 1 end
function Vlasov:isPositive(q)
   if q[0] > 0.0 then
      return true
   else
      return false
   end
end

-- flux in direction dir
function Vlasov:flux(dir, qIn, fOut)
   assert(false, "Vlasov:flux: NYI!")
end

-- Riemann problem for Vlasov equation
function Vlasov:rp(dir, delta, ql, qr, waves, s)
   assert(false, "Vlasov:rp: NYI!")
end

-- Compute q-fluctuations
function Vlasov:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "Vlasov:qFluctuations: NYI!")
end

-- Maximum wave speed
function Vlasov:maxSpeed(dir, w, dx, q)
   return 0.0
end       
-- Volume integral term for use in DG scheme
function Vlasov:volTerm(w, dx, idx, q, out)
   -- streaming term
   local cflFreqStream = self._volStreamUpdate(w:data(), dx:data(), q:data(), out:data())
   local cflFreqForce = 0.0
   -- force term
   if self._hasForceTerm then
      -- set pointer to EM field and scale by q/m
      self._emField:fill(self._emIdxr(idx), self._emPtr) -- get pointer to EM field
      rescaleEmField(self._qbym, self._emPtr, self._emAccel) -- multiply EM field by q/m
      -- compute vol term
      cflFreqForce = self._volForceUpdate(w:data(), dx:data(), self._emAccel:data(), q:data(), out:data())
   end
   return cflFreqStream+cflFreqForce
end

-- Surface integral term for use in DG scheme
function Vlasov:surfTerm(dir, w, dx, maxs, idxl, idxr, ql, qr, outl, outr)
   local amax = 0.0
   if dir <= self._cdim then
      -- streaming term (note that surface streaming kernels don't
      -- return max speed)
      self._surfStreamUpdate[dir](
	 w:data(), dx:data(), ql:data(), qr:data(), outl:data(), outr:data())
   else
      if self._hasForceTerm then
	 -- force term
	 -- set pointer to EM field and scale by q/m
	 self._emField:fill(self._emIdxr(idxl), self._emPtr) -- get pointer to EM field
	 rescaleEmField(self._qbym, self._emPtr, self._emAccel) -- multiply EM field by q/m
	 amax = self._surfForceUpdate[dir-self._cdim](
	    w:data(), dx:data(), maxs, self._emAccel:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return amax
end

function Vlasov:setAuxFields(auxFields)
   if self._hasForceTerm then -- (no fields for neutral particles)
      -- single aux field that has the full EM field
      self._emField = auxFields[1]

      if self._isFirst then
	 -- allocate pointers to field object
	 self._emPtr = self._emField:get(1)
	 self._emIdxr = self._emField:genIndexer()
	 self._emAccel = Lin.Vec(self._emField:numComponents())
	 self._isFirst = false -- no longer first time
      end
   end
end

function Vlasov:pdim() return self._pdim end
function Vlasov:cdim() return self._cdim end
function Vlasov:vdim() return self._vdim end
function Vlasov:qbym() return self._qbym end
function Vlasov:hasForceTerm() return self._hasForceTerm end

function Vlasov:volStreamUpdate()
   return self._volStreamUpdate
end
function Vlasov:surfStreamUpdate()
   return self._surfStreamUpdate
end

function Vlasov:volForceUpdate()
   return self._volForceUpdate
end
function Vlasov:surfForceUpdate()
   return self._surfForceUpdate
end

return Vlasov
