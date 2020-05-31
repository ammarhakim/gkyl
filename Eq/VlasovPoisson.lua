-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov-Poisson equation on a rectangular mesh.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local EqBase               = require "Eq.EqBase"
local Lin                  = require "Lib.Linalg"
local Proto                = require "Lib.Proto"
local VlasovPoissonModDecl = require "Eq.vlasovPoissonData.VlasovPoissonModDecl"
local ffi                  = require "ffi"
local ffiC                 = ffi.C
local xsys                 = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

ffi.cdef [[ 
  typedef struct VlasovPoisson VlasovPoisson;
  Vlasov* new_VlasovPoisson(unsigned cdim, unsigned vdim, unsigned polyOrder, const char* basisType, double qbym, bool hasForceTerm);
  void setAuxFields(GkylEquation_t *eq, GkylCartField_t *phi);
]]

-- Vlasov-Poisson equation on a rectangular mesh.
local VlasovPoisson = Proto(EqBase)

-- ctor
function VlasovPoisson:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.VlasovPoisson: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Eq.VlasovPoisson: Must specify configuration-space basis functions to use using 'confBasis'")
   
   local charge = assert(tbl.charge, "Eq.VlasovPoisson: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "Eq.VlasovPoisson: must specify mass using 'mass' ")

   self._qbym = charge/mass   -- Only the q/m ratio is ever needed.

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- Functions to perform streaming updates.
   self._volUpdate        = VlasovPoissonModDecl.selectVolStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovPoissonModDecl.selectSurfStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   -- Check if we have a electric field.
   local hasElcField = xsys.pickBool(tbl.hasElectricField, true)

   self._hasForceTerm = false   -- Flag to indicate if we have any force terms at all.
   if hasElcField then self._hasForceTerm = true end

   self._volForceUpdate, self._surfForceUpdate = nil, nil
   if self._hasForceTerm then
      -- Functions to perform force updates.
      self._volUpdate = VlasovPoissonModDecl.selectVolElc(
         self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfForceUpdate = VlasovPoissonModDecl.selectSurfElc(
         self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end

   -- Electrostatic potential (phi) object and pointers to cell values.
   self._phi = nil
   self._phiPtr, self._phiIdxr = nil, nil

   -- Flag to indicate if we are being called for first time.
   self._isFirst = true

   if GKYL_HAVE_CUDA then
      self:initDevice()
   end
end

function VlasovPoisson:initDevice()
   local eq = ffiC.new_VlasovPoisson(self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._phaseBasis:id(), self._qbym, self._hasForceTerm) 
   local sz = sizeof(eq)
   self._onDevice, err = cuda.Malloc(sz)
   cuda.Memcpy(self._onDevice, eq, sz, cuda.MemcpyHostToDevice)
end

function VlasovPoisson:numEquations() return 1 end
function VlasovPoisson:numWaves() return 1 end
function VlasovPoisson:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- Flux in direction dir.
function VlasovPoisson:flux(dir, qIn, fOut)
   assert(false, "VlasovPoisson:flux: NYI!")
end

-- Riemann problem for VlasovPoisson equation.
function VlasovPoisson:rp(dir, delta, ql, qr, waves, s)
   assert(false, "VlasovPoisson:rp: NYI!")
end

-- Compute q-fluctuations.
function VlasovPoisson:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "VlasovPoisson:qFluctuations: NYI!")
end

-- Maximum wave speed.
function VlasovPoisson:maxSpeed(dir, w, dx, q)
   assert(false, "VlasovPoisson:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function VlasovPoisson:volTerm(w, dx, idx, q, out)
   -- Volume term if has force.
   local cflFreq = 0.0
   if self._hasForceTerm then
      self._phi:fill(self._phiIdxr(idx), self._phiPtr)   -- Get pointer to electrostatic potential.
      cflFreq = self._volUpdate(w:data(), dx:data(), self._phiPtr:data(), q:data(), out:data())
   else
      -- If no force, only update streaming term.
      cflFreq = self._volUpdate(w:data(), dx:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function VlasovPoisson:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local amax = 0.0
   if dir <= self._cdim then
      -- Streaming term (note that surface streaming kernels don't return max speed).
      self._surfStreamUpdate[dir](
	 wl:data(), wr:data(), dxl:data(), dxr:data(), ql:data(), qr:data(), outl:data(), outr:data())
   else
      if self._hasForceTerm then
	 -- Force term.
	 self.phiField:fill(self._phiIdxr(idxl), self._phiPtr)   -- Get pointer to electrostatic potential.
	 amax = self._surfForceUpdate[dir-self._cdim](wl:data(), wr:data(), dxl:data(), dxr:data(), maxs,
	    self._phiPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   end
   return amax
end

function VlasovPoisson:setAuxFields(auxFields)
   if self._hasForceTerm then

      self._phi = auxFields[1]

      if self._isFirst then
	 self._phiPtr  = self._phiField:get(1)
	 self._phiIdxr = self._phiField:genIndexer()
	 self._isFirst = false
      end
   end
end

function VlasovPoisson:setAuxFieldsOnDevice(auxFields)
   if self._hasForceTerm then -- (no fields for neutral particles)
      self._phi = auxFields[1]
      ffiC.setAuxFields(self._onDevice, self._phi._onDevice)
   end
end

function VlasovPoisson:pdim() return self._pdim end
function VlasovPoisson:cdim() return self._cdim end
function VlasovPoisson:vdim() return self._vdim end
function VlasovPoisson:qbym() return self._qbym end
function VlasovPoisson:hasForceTerm() return self._hasForceTerm end

return VlasovPoisson
