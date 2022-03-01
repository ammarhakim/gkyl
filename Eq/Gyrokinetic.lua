-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic equation using Hamiltonian formulation.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local EqBase     = require "Eq.EqBase"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local xsys       = require "xsys"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local ffi           = require "ffi"
local ffiC          = ffi.C

ffi.cdef [[ 

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_gyrokinetic_auxfields {
  const struct gkyl_array *bmag; // Pointer to magnetic field magnitude.
  const struct gkyl_array *jacobtot_inv; // Pointer to 1/(conf Jacobian * gyro-center coords Jacobian).
  const struct gkyl_array *cmag; // Pointer to parallel gradient coefficient.
  const struct gkyl_array *b_i; // Pointer to covariant components of B-field unit vector.
  const struct gkyl_array *phi; // Pointer to electrostatic potential.
  const struct gkyl_array *apar; // Pointer to A_\parallel.
  const struct gkyl_array *apardot; // Pointer to d(A_parallel)/dt.
};

/**                                                                                         
 * Create a new Gyrokinetic equation object.                                                     
 *                                                                                          
 * @param cbasis Configuration space basis functions                                        
 * @param pbasis Phase-space basis functions                                                
 * @param conf_range Configuration space range for use in indexing EM field                 
 * @return Pointer to Gyrokinetic equation object                                                
 */                                                                                         
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_new(const struct gkyl_basis* cbasis,                     
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const double charge, const double mass, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. geometry & EM fields) needed in computing
 * gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin);

]]

local Gyrokinetic = Proto(EqBase)

-- ctor
function Gyrokinetic:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "Gyrokinetic: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "Gyrokinetic: must specify a phaseBasis")
   self._confGrid  = assert(tbl.confGrid, "Gyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "Gyrokinetic: must specify confBasis")

   self._confRange = tbl.confRange
   if self._confRange then assert(self._confRange:isSubRange()==1, "Eq.Gyrokinetic: confRange must be a sub-range") end

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)
   self._positivity        = xsys.pickBool(tbl.positivity,false)

   self.Bvars = tbl.Bvars

   self.geoType = tbl.geometry and tbl.geometry or "SimpleHelical"

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()

   self._zero = ffiC.gkyl_dg_gyrokinetic_new(self._confBasis._zero, self._basis._zero, self._confRange, self.charge, self.mass, GKYL_USE_GPU or 0)

   self._volTerm  = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars, self.geoType)
   self._surfTerm = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self._positivity, self.Bvars, self.geoType)

   -- Select the appropriate volume and surface term functions to call.
   self.volTermFunc  = function(w, dx, idx, f, out) return Gyrokinetic["volTerm" .. self.geoType](self, w, dx, idx, f, out) end
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return Gyrokinetic["surfTerm" .. self.geoType](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end

   -- For sheath BCs.
   if tbl.hasSheathBCs then
      self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(nm, self._cdim, self._vdim, p)
   end

   if self._isElectromagnetic then
      self.emMod = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1}
      }
      self.emModPtrL = self.emMod:get(1)
      self.emModPtrR = self.emMod:get(1)
      self.emModIdxr = self.emMod:genIndexer()
      self.emMod:clear(0.0)
   end

   if self._zero then self._auxfieldsC = ffi.new("struct gkyl_dg_gyrokinetic_auxfields") end

   -- For gyroaveraging.
   self.gyavgSlvr = tbl.gyavgSlvr
   if self.gyavgSlvr then
      self._gyavg = true
      self.phiGy = {}
      local nmu = self._grid:numCells(self._ndim)
      for i = 1, nmu do
         self.phiGy[i] = DataStruct.Field {
            onGrid        = self._confGrid,
            numComponents = self._confBasis:numBasis(),
            ghost         = {1, 1}
         }
      end
   end

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function Gyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   -- Get electromagnetic terms.
   self.apar = potentials.apar
   self.dApardt = potentials.dApardt
   self.dApardtProv = auxFields[3]

   if self._zero then
      self._auxfieldsC.bmag = geo.bmag._zero
      self._auxfieldsC.jacobtot_inv = geo.jacobTotInv._zero
      self._auxfieldsC.cmag = geo.cmag._zero
      self._auxfieldsC.b_i = geo.b_i._zero
      self._auxfieldsC.phi = self.phi._zero
      self._auxfieldsC.apar = self.apar._zero
      self._auxfieldsC.apardot = self.dApardt._zero
      ffiC.gkyl_gyrokinetic_set_auxfields(self._zero, self._auxfieldsC)
      return
   end

   -- Get magnetic geometry fields.
   self.bmag = geo.bmag
   self.cmag = geo.cmag
   if self.geoType == "SimpleHelical" then
      self.bmagInv = geo.bmagInv
      self.bdriftX = geo.bdriftX
      self.bdriftY = geo.bdriftY
   elseif self.geoType == "GenGeo" then
      self.b_x = geo.b_x
      self.b_y = geo.b_y
      self.b_z = geo.b_z
      self.jacobTotInv = geo.jacobTotInv
   end
   self.phiWall = geo.phiWall  -- For sheath BCs.



   if self._isFirst then
      -- Allocate pointers and indexers to field objects.

      -- Potentials.
      self.phiPtr  = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      if self._isElectromagnetic then
         self.aparPtr        = self.apar:get(1)
         self.aparLPtr       = self.apar:get(1)
         self.dApardtPtr     = self.dApardt:get(1)
         self.dApardtProvPtr = self.dApardtProv:get(1)
         self.aparIdxr       = self.apar:genIndexer()
         self.dApardtIdxr    = self.dApardt:genIndexer()
      end

      -- For gyroaveraging.
      if self._gyavg then
         self.phiGyPtr = {}
         self.phiGyIdxr = {}
         for i=1,self._grid:numCells(self._ndim) do
            self.phiGyPtr[i]  = self.phiGy[i]:get(1)
            self.phiGyIdxr[i] = self.phiGy[i]:genIndexer()
         end
      end

      -- Geometry.
      self.bmagPtr     = self.bmag:get(1)
      self.cmagPtr     = self.cmag:get(1)
      self.phiWallPtr  = self.phiWall:get(1)
      self.bmagIdxr    = self.bmag:genIndexer()
      self.cmagIdxr    = self.cmag:genIndexer()
      self.phiWallIdxr = self.phiWall:genIndexer()
      if self.geoType == "SimpleHelical" then
         self.bmagInvPtr  = self.bmagInv:get(1)
         self.bdriftXPtr  = self.bdriftX:get(1)
         self.bdriftYPtr  = self.bdriftY:get(1)
         self.bmagInvIdxr = self.bmagInv:genIndexer()
         self.bdriftXIdxr = self.bdriftX:genIndexer()
         self.bdriftYIdxr = self.bdriftY:genIndexer()
      elseif self.geoType == "GenGeo" then
         self.jacobTotInvPtr  = self.jacobTotInv:get(1)
         self.b_xPtr         = self.b_x:get(1)
         self.b_yPtr         = self.b_y:get(1)
         self.b_zPtr         = self.b_z:get(1)
         self.jacobTotInvIdxr = self.jacobTotInv:genIndexer()
         self.b_xIdxr        = self.b_x:genIndexer()
         self.b_yIdxr        = self.b_y:genIndexer()
         self.b_zIdxr        = self.b_z:genIndexer()
      end

      self._isFirst = false -- No longer first time.
   end
end

function Gyrokinetic:setAuxFieldsOnDevice(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   -- Get electromagnetic terms.
   self.apar = potentials.apar
   self.dApardt = potentials.dApardt
   if self._zero then
      self._auxfieldsC.bmag = geo.bmag._zeroDot
      self._auxfieldsC.jacobtot_inv = geo.jacobTotInv._zeroDot
      self._auxfieldsC.cmag = geo.cmag._zeroDot
      self._auxfieldsC.b_i = geo.b_i._zeroDot
      self._auxfieldsC.phi = self.phi._zeroDot
      self._auxfieldsC.apar = self.apar._zeroDot
      self._auxfieldsC.apardot = self.dApardt._zeroDot
      ffiC.gkyl_gyrokinetic_set_auxfields(self._zero, self._auxfieldsC)
      return
   end
end

-- Volume integral term for use in DG scheme.
function Gyrokinetic:volTerm(w, dx, idx, f, out)
   return self.volTermFunc(w, dx, idx, f, out)
end
function Gyrokinetic:volTermSimpleHelical(w, dx, idx, f, out)
   local tmStart = Time.clock()
   if self._gyavg then
      local idmu = idx[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idx), self.phiGyPtr[idmu])
      self.phiPtr = self.phiGyPtr[idmu]
   else
      self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   end
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idx), self.bmagInvPtr)
   self.cmag:fill(self.cmagIdxr(idx), self.cmagPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idx), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idx), self.bdriftYPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idx), self.aparPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idx), self.dApardtProvPtr)
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.cmagPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtProvPtr:data(), f:data(), out:data())
   else
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.cmagPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), f:data(), out:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end
function Gyrokinetic:volTermGenGeo(w, dx, idx, f, out)
   local tmStart = Time.clock()
   if self._gyavg then 
      local idmu = idx[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idx), self.phiGyPtr[idmu])
      self.phiPtr = self.phiGyPtr[idmu]
   else
      self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   end
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idx), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idx), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idx), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idx), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idx), self.b_zPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idx), self.aparPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idx), self.dApardtProvPtr)
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtProvPtr:data(), f:data(), out:data())
   else 
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), f:data(), out:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function Gyrokinetic:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function Gyrokinetic:surfTermSimpleHelical(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   if self._gyavg then
      local idmu = idxr[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idxr), self.phiGyPtr[idmu])
      self.phiPtr = self.phiGyPtr[idmu]
   else
      self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   end
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftYPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)
     self.emMod:fill(self.emModIdxr(idxl), self.emModPtrL)
     self.emMod:fill(self.emModIdxr(idxr), self.emModPtrR)
     res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.cmagPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), self.emModPtrL:data(), self.emModPtrR:data())
   else
     res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.cmagPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end
function Gyrokinetic:surfTermGenGeo(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   if self._gyavg then 
      local idmu = idxr[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idxr), self.phiGyPtr[idmu])
      self.phiPtr = self.phiGyPtr[idmu]
   else
      self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   end
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
     self.apar:fill(self.aparIdxr(idxl), self.aparLPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)
     self.emMod:fill(self.emModIdxr(idxl), self.emModPtrL)
     self.emMod:fill(self.emModIdxr(idxr), self.emModPtrR)
     res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.aparLPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), self.emModPtrL:data(), self.emModPtrR:data())
   else 
     res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function Gyrokinetic:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f, fRefl)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.phiWall:fill(self.phiWallIdxr(idx), self.phiWallPtr)
   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, 
                                     self.phiPtr:data(), self.phiWallPtr:data(), f:data(), fRefl:data())
end

local GyrokineticStep2 = Proto(EqBase)
-- ctor
function GyrokineticStep2:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "GyrokineticStep2: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "GyrokineticStep2: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "GyrokineticStep2: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "GyrokineticStep2: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "GyrokineticStep2: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   self._positivity = xsys.pickBool(tbl.positivity,false)

   self.Bvars = tbl.Bvars

   self.geoType = tbl.geometry and tbl.geometry or "SimpleHelical"

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = GyrokineticModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p, self.geoType)
   self._surfTerm = GyrokineticModDecl.selectStep2Surf(nm, self._cdim, self._vdim, p, self._positivity, self.Bvars, self.geoType)

   -- Select the appropriate volume and surface term functions to call.
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return GyrokineticStep2["surfTerm" .. self.geoType](self,dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end

   self._isFirst = true
end

function GyrokineticStep2:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.
   --local potentialsProv = auxFields[3]

   -- Get phi, Apar, and dApar/dt.
   self.phi         = potentials.phi
   self.apar        = potentials.apar
   self.dApardt     = potentials.dApardt
   self.dApardtProv = auxFields[3]

   -- Get magnetic geometry fields.
   self.bmag = geo.bmag
   self.cmag = geo.cmag
   if self.geoType == "SimpleHelical" then
      self.bmagInv = geo.bmagInv
      self.bdriftX = geo.bdriftX
      self.bdriftY = geo.bdriftY
   elseif self.geoType == "GenGeo" then
      self.b_x = geo.b_x
      self.b_y = geo.b_y
      self.b_z = geo.b_z
      self.jacobTotInv = geo.jacobTotInv
   end

   if self._isFirst then
      -- Allocate pointers and indexers to field objects.

      -- Potentials.
      self.phiPtr         = self.phi:get(1)
      self.phiIdxr        = self.phi:genIndexer()
      self.aparPtr        = self.apar:get(1)
      self.aparLPtr       = self.apar:get(1)
      self.dApardtPtr     = self.dApardt:get(1)
      self.dApardtProvPtr = self.dApardtProv:get(1)
      self.aparIdxr       = self.apar:genIndexer()
      self.dApardtIdxr    = self.dApardt:genIndexer()

      -- Geometry.
      self.bmagPtr  = self.bmag:get(1)
      self.cmagPtr  = self.cmag:get(1)
      self.bmagIdxr = self.bmag:genIndexer()
      self.cmagIdxr = self.cmag:genIndexer()
      if self.geoType == "SimpleHelical" then
         self.bmagInvPtr  = self.bmagInv:get(1)
         self.bdriftXPtr  = self.bdriftX:get(1)
         self.bdriftYPtr  = self.bdriftY:get(1)
         self.bmagInvIdxr = self.bmagInv:genIndexer()
         self.bdriftXIdxr = self.bdriftX:genIndexer()
         self.bdriftYIdxr = self.bdriftY:genIndexer()
      elseif self.geoType == "GenGeo" then
         self.jacobTotInvPtr  = self.jacobTotInv:get(1)
         self.b_xPtr         = self.b_x:get(1)
         self.b_yPtr         = self.b_y:get(1)
         self.b_zPtr         = self.b_z:get(1)
         self.jacobTotInvIdxr = self.jacobTotInv:genIndexer()
         self.b_xIdxr        = self.b_x:genIndexer()
         self.b_yIdxr        = self.b_y:genIndexer()
         self.b_zIdxr        = self.b_z:genIndexer()
      end

      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function GyrokineticStep2:volTerm(w, dx, idx, f, out)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.dApardt:fill(self.dApardtIdxr(idx), self.dApardtPtr)
   return self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.dApardtPtr:data(), f:data(), out:data())
end

-- Surface integral term for use in DG scheme 
-- NOTE: only vpar direction for this term
function GyrokineticStep2:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function GyrokineticStep2:surfTermSimpleHelical(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftYPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)

   local res = self._surfTerm(self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.cmagPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), nil, nil)

   return res
end
function GyrokineticStep2:surfTermGenGeo(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.apar:fill(self.aparIdxr(idxl), self.aparLPtr)
   self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)

   local res = self._surfTerm(self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.aparLPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), nil, nil)

   return res
end

return {GkEq = Gyrokinetic, GkEqStep2 = GyrokineticStep2} 
