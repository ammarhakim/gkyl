-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local EqBase = require "Eq.EqBase"
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
local Vlasov = Proto(EqBase)

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
   self._volStreamUpdate = VlasovModDecl.selectVolStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovModDecl.selectSurfStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   -- check if we have a electric and magnetic field
   local hasElcField = xsys.pickBool(tbl.hasElectricField, true)
   local hasMagField = xsys.pickBool(tbl.hasMagneticField, true)

   self._hasForceTerm = false -- flag to indicate if we have any force terms at all
   if hasElcField or hasMagField then
      self._hasForceTerm = true
   end

   self._volForceUpdate, self._surfForceUpdate = nil, nil
   if self._hasForceTerm then
      -- functions to perform force updates
      if hasMagField then 
	 self._volUpdate = VlasovModDecl.selectVolElcMag(
	    self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
	 self._surfForceUpdate = VlasovModDecl.selectSurfElcMag(
	    self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      else
	 assert(false, "Vlasov: Pure ES kernels NYI!")
      end
   end

   -- EM field object and pointers to cell values
   self._emField = nil
   -- (these will be set on the first call to setAuxFields() method)
   self._emPtr, self._emIdxr = nil, nil
   self._emAccel = nil

   self._hasConstGrav = false
   self._gravAccelIdx, self._gravAccelCoeff = 1, 0.0
   -- check if we have gravity terms
   if tbl.constGravity then
      self._hasConstGrav = true
      -- gravAccelIdx is the index location to add gravitational
      -- acceleration term
      self._gravAccelIdx = 1+self._confBasis:numBasis()*(tbl.constGravity.dir-1)
      self._gravAccelCoeff = math.sqrt(2^self._cdim)*tbl.constGravity.accel
   end

   -- flag to indicate if we are being called for first time
   self._isFirst = true
end

-- Methods
function Vlasov:numEquations() return 1 end
function Vlasov:numWaves() return 1 end
function Vlasov:isPositive(q)
   if q[1] > 0.0 then
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
   assert(false, "Vlasov:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme
function Vlasov:volTerm(w, dx, idx, q, out)
   -- volume term if has force
   local cflFreq = 0.0
   if self._hasForceTerm then
      self._emField:fill(self._emIdxr(idx), self._emPtr) -- get pointer to EM field
      rescaleEmField(self._qbym, self._emPtr, self._emAccel) -- multiply EM field by q/m
      if self._hasConstGrav then
	 self._emAccel[self._gravAccelIdx] = self._emAccel[self._gravAccelIdx]+self._gravAccelCoeff
      end
      cflFreq = self._volUpdate(w:data(), dx:data(), self._emAccel:data(), q:data(), out:data())
   else
      -- if no force, only update streaming term
      cflFreq = self._volUpdate(w:data(), dx:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme
function Vlasov:surfTerm(dir, dt, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local amax = 0.0
   if dir <= self._cdim then
      -- streaming term (note that surface streaming kernels don't
      -- return max speed)
      self._surfStreamUpdate[dir](
	 wl:data(), wr:data(), dxl:data(), dxr:data(), ql:data(), qr:data(), outl:data(), outr:data())
   else
      if self._hasForceTerm then
	 -- force term
	 self._emField:fill(self._emIdxr(idxl), self._emPtr) -- get pointer to EM field
	 rescaleEmField(self._qbym, self._emPtr, self._emAccel) -- multiply EM field by q/m
	 if self._hasConstGrav then
	    self._emAccel[self._gravAccelIdx] = self._emAccel[self._gravAccelIdx]+self._gravAccelCoeff
	 end	 
	 amax = self._surfForceUpdate[dir-self._cdim](
	    wl:data(), wr:data(), dxl:data(), dxr:data(), maxs,
	    self._emAccel:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return amax
end

function Vlasov:setAuxFields(auxFields)
   if self._hasForceTerm then -- (no fields for neutral particles)
      -- single aux field that has the full EM field
      self._emField = auxFields[1]

      if self._isFirst then
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

return Vlasov
