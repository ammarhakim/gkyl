-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov equation on a rectangular mesh.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local CartField     = require "DataStruct.CartField"
local EqBase        = require "Eq.EqBase"
local Lin           = require "Lib.Linalg"
local Proto         = require "Lib.Proto"
local VlasovModDecl = require "Eq.vlasovData.VlasovModDecl"
local ffi           = require "ffi"
local ffiC          = ffi.C
local xsys          = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local cuda = nil
if GKYL_HAVE_CUDA then cuda = require "Cuda.RunTime" end

ffi.cdef [[ 
// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_PHI, // Poisson (only phi)  
  GKYL_FIELD_PHI_A, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL // no field is present
};

/**                                                                                         
 * Create a new Vlasov equation object.                                                     
 *                                                                                          
 * @param cbasis Configuration space basis functions                                        
 * @param pbasis Phase-space basis functions                                                
 * @param conf_range Configuration space range for use in indexing EM field                 
 * @param field_id enum to determine what type of EM fields (Vlasov-Maxwell vs. neutrals)   
 * @return Pointer to Vlasov equation object                                                
 */                                                                                         
struct gkyl_dg_eqn* gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis,                     
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, enum gkyl_field_id field_id, bool use_gpu);
  
/**
 * Set the q/m*EM field needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param qmem Pointer to EM field scaled by q/m,
 */
void gkyl_vlasov_set_qmem(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem);
]]

-- Vlasov equation on a rectangular mesh.
local Vlasov = Proto(EqBase)

-- ctor
function Vlasov:init(tbl)

   self._phaseBasis = assert(
      tbl.phaseBasis, "Eq.Vlasov: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Eq.Vlasov: Must specify configuration-space basis functions to use using 'confBasis'")

   self._confRange = tbl.confRange
   if self._confRange then assert(self._confRange:isSubRange()==1, "Eq.Vlasov: confRange must be a sub-range") end
   
   local charge = assert(tbl.charge, "Eq.Vlasov: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Eq.Vlasov: must specify mass using 'mass' ")

   self._qbym = charge/mass   -- Only q/m ratio is ever needed.

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- Functions to perform streaming updates.
   self._volUpdate = VlasovModDecl.selectVolStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovModDecl.selectSurfStream(
      self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   -- Check if we have an electric and magnetic field.
   local hasElcField    = xsys.pickBool(tbl.hasElectricField, true)
   local hasMagField    = xsys.pickBool(tbl.hasMagneticField, true)
   local hasExtForce    = xsys.pickBool(tbl.hasExtForce, true)
   self._plasmaMagField = xsys.pickBool(tbl.plasmaMagField, true)

   if hasElcField and hasMagField then 
      self._fieldId = 0 
   elseif hasElcField then
      self._fieldId = 1
   else
      self._fieldId = 3
   end
   self._zero = ffiC.gkyl_dg_vlasov_new(self._confBasis._zero, self._phaseBasis._zero, self._confRange, self._fieldId, GKYL_USE_GPU or 0)

   -- Option to perform only the force volume update (e.g. for stochastic forces).
   local onlyForceUpdate = xsys.pickBool(tbl.onlyForceUpdate, false)

   self._hasForceTerm = false   -- Flag to indicate if we have any force terms at all.
   -- Turn on force term if there's an electric or magnetic field or if there's a body force.
   if hasElcField or hasMagField or self._plasmaMagField or hasExtForce then
      self._hasForceTerm = true
   end

   self._onlyForceUpdate = false   -- Flag to indicate if updating force separately from streaming.
   if onlyForceUpdate then self._onlyForceUpdate = true end

   self._surfForceUpdate = nil
   -- numVelFlux used for selecting which type of numerical flux function to use for velocity space update
   -- default is "penalty," supported options: "penalty," "recovery"
   self._numVelFlux = tbl.numVelFlux and tbl.numVelFlux or "penalty"
   if self._hasForceTerm then
      if self._qbym == 0. and hasExtForce then
         self._volUpdate = VlasovModDecl.selectVolNeutral(
            self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
         self._surfForceUpdate = VlasovModDecl.selectSurfNeutral(
            self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
      else
         -- Functions to perform force updates.
         if self._plasmaMagField then 
            self._volUpdate = VlasovModDecl.selectVolElcMag(
               self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
            self._surfForceUpdate = VlasovModDecl.selectSurfElcMag(
               self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._numVelFlux)
         elseif self._onlyForceUpdate then
   	 self._volUpdate = VlasovModDecl.selectVolForce(
   	    self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
            self._surfForceUpdate = VlasovModDecl.selectSurfElcMag(
               self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder(), self._numVelFlux)
         else
            self._volUpdate = VlasovModDecl.selectVolPhi(hasMagField,
               self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
            self._surfForceUpdate = VlasovModDecl.selectSurfPhi(hasMagField,
               self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
         end
      end
   end

   -- EM field object and pointers to cell values.
   self._emField = nil
   -- (these will be set on the first call to setAuxFields() method).
   self._emPtr, self._emIdxr = nil, nil

   -- Electrostatic potential (for Vlasov-Poisson).
   self._phiField = nil
   self._phiPtr, self._phiIdxr = nil, nil

   -- Flag to indicate if we are being called for first time.
   self._isFirst = true
end

-- Methods.
function Vlasov:numEquations() return 1 end
function Vlasov:numWaves() return 1 end
function Vlasov:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- Flux in direction dir.
function Vlasov:flux(dir, qIn, fOut)
   assert(false, "Vlasov:flux: NYI!")
end

-- Riemann problem for Vlasov equation.
function Vlasov:rp(dir, delta, ql, qr, waves, s)
   assert(false, "Vlasov:rp: NYI!")
end

-- Compute q-fluctuations.
function Vlasov:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "Vlasov:qFluctuations: NYI!")
end

-- Maximum wave speed.
function Vlasov:maxSpeed(dir, w, dx, q)
   assert(false, "Vlasov:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function Vlasov:volTerm(w, dx, idx, q, out)
   -- Volume term if has force.
   local cflFreq = 0.0
   if self._hasForceTerm then
      self._emField:fill(self._emIdxr(idx), self._emPtr)   -- Get pointer to EM field.
      if self._plasmaMagField then 
         cflFreq = self._volUpdate(w:data(), dx:data(), self._emPtr:data(), q:data(), out:data())
      else
         self._phiField:fill(self._phiIdxr(idx), self._phiPtr)   -- Get pointer to the electrostatic potential.
         cflFreq = self._volUpdate(w:data(), dx:data(), self._qbym, self._phiPtr:data(), self._emPtr:data(), q:data(), out:data())
      end
   else
      -- If no force, only update streaming term.
      cflFreq = self._volUpdate(w:data(), dx:data(), q:data(), out:data())
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function Vlasov:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   local amax = 0.0
   if dir <= self._cdim then
      if not self._onlyForceUpdate then
         -- Streaming term (note that surface streaming kernels don't return max speed).
         self._surfStreamUpdate[dir](
           wl:data(), wr:data(), dxl:data(), dxr:data(), ql:data(), qr:data(), outl:data(), outr:data())
      end
   else
      if self._hasForceTerm then
	 -- Force term.
	 self._emField:fill(self._emIdxr(idxl), self._emPtr)   -- Get pointer to EM field.
         if self._plasmaMagField then 
	    amax = self._surfForceUpdate[dir-self._cdim](
	       wl:data(), wr:data(), dxl:data(), dxr:data(), maxs,
	       self._emPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         else
	   self._phiField:fill(self._phiIdxr(idxl), self._phiPtr)   -- Get pointer to the electrostatic potential.
	    amax = self._surfForceUpdate[dir-self._cdim](
	       wl:data(), wr:data(), dxl:data(), dxr:data(), maxs, self._qbym, self._phiPtr:data(),
	       self._emPtr:data(), ql:data(), qr:data(), outl:data(), outr:data())
         end
      end
   end
   return amax
end

function Vlasov:setAuxFields(auxFields)
   if self._hasForceTerm then   -- No fields for neutral particles.
      self._emField = auxFields[1]   -- Single aux field that has the full EM field
      if not self._plasmaMagField then self._phiField = auxFields[2] end   -- Electrostatic potential.

      if self._isFirst then
         self._emPtr  = self._emField:get(1)
         self._emIdxr = self._emField:genIndexer()
         if not self._plasmaMagField then 
            self._phiPtr  = self._phiField:get(1)
            self._phiIdxr = self._phiField:genIndexer()
         end
         self._isFirst = false   -- No longer first time.
      end
      if self._zero then
         ffiC.gkyl_vlasov_set_qmem(self._zero, self._emField._zero)
      end
   end
end

function Vlasov:setAuxFieldsOnDevice(auxFields)
   if self._hasForceTerm then   -- No fields for neutral particles.
      -- Single aux field that has the full EM field.
      self._emField = auxFields[1]
      ffiC.gkyl_vlasov_set_qmem(self._zero, self._emField._zeroDevice);
   end
end

function Vlasov:pdim() return self._pdim end
function Vlasov:cdim() return self._cdim end
function Vlasov:vdim() return self._vdim end
function Vlasov:qbym() return self._qbym end
function Vlasov:hasForceTerm() return self._hasForceTerm end

return Vlasov
