-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local BoundaryCondition = require "Updater.BoundaryCondition"
local Proto = require "Lib.Proto"
local VlasovModDecl = require "Updater.vlasovData.VlasovModDecl"
local ffi = require "ffi"
local xsys = require "xsys"

-- Vlasov equation on a rectangular mesh
local VlasovRect = Proto()

-- ctor
function VlasovRect:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.VlasovRect: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Eq.VlasovRect: Must specify configuration-space basis functions to use using 'confBasis'")
   
   local charge = assert(tbl.charge, "Eq.VlasovRect: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Eq.VlasovRect: must specify mass using 'mass' ")

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
      self._volForceUpdate = VlasovModDecl.selectVolElcMag(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfForceUpdate = VlasovModDecl.selectSurfElcMag(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   else
      self._volForceUpdate = VlasovModDecl.selectVolElc(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      self._surfForceUpdate = VlasovModDecl.selectSurfElc(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   end   
end

-- Methods
function VlasovRect:numEquations() return 1 end
function VlasovRect:numWaves() return 1 end
function VlasovRect:isPositive(q)
   if q[0] > 0.0 then
      return true
   else
      return false
   end
end

-- flux in direction dir
function VlasovRect:flux(dir, qIn, fOut)
   assert(false, "VlasovRect:flux: NYI!")
end

-- Riemann problem for Vlasov equation
function VlasovRect:rp(dir, delta, ql, qr, waves, s)
   assert(false, "VlasovRect:rp: NYI!")
end

-- Compute q-fluctuations
function VlasovRect:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "VlasovRect:qFluctuations: NYI!")
end

-- Maximum wave speed
function VlasovRect:maxSpeed(dir, w, dx, q)
   return 0.0
end       
-- Volume integral term for use in DG scheme
function VlasovRect:volTerm(w, dx, idx, q, out)
   -- TODO
end
-- Surface integral term for use in DG scheme
function VlasovRect:surfTerm(dir, w, dx, maxs, idxl, idxr, ql, qr, outl, outr)
   -- TODO
end

function VlasovRect:setAuxFields(auxFields)
end

return VlasovRect
