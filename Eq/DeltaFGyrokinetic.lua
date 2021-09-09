-- Gkyl ------------------------------------------------------------------------
--
-- DeltaFGyrokinetic equation using Hamiltonian formulation.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local EqBase     = require "Eq.EqBase"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local xsys       = require "xsys"
local DeltaFGyrokineticModDecl = require "Eq.gkData.DeltaFGyrokineticModDecl"

local DeltaFGyrokinetic = Proto(EqBase)

-- ctor
function DeltaFGyrokinetic:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "DeltaFGyrokinetic: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "DeltaFGyrokinetic: must specify a phaseBasis")
   self._confGrid  = assert(tbl.confGrid, "DeltaFGyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "DeltaFGyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "DeltaFGyrokinetic: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "DeltaFGyrokinetic: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   assert(tbl.hasPhi==true, "DeltaFGyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)
   self._positivity        = xsys.pickBool(tbl.positivity,false)

   --self.f0 = assert(tbl.f0, "DeltaFGyrokinetic: must specify a background distribution function using 'f0' ")
   self.linear = tbl.linear

   self.Bvars = tbl.Bvars

   self.geoType = "GenGeo"

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = DeltaFGyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars, self.linear)
   self._surfTerm = DeltaFGyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars, self.linear)

   -- Select the appropriate volume and surface term functions to call.
   self.volTermFunc  = function(w, dx, idx, f, out) return DeltaFGyrokinetic["volTerm" .. self.geoType](self, w, dx, idx, f, out) end
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return DeltaFGyrokinetic["surfTerm" .. self.geoType](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
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

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function DeltaFGyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.
   self.f0 = auxFields[4]

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi
   if self._gyavg then 
      self.gyavgSlvr:advance(0, {self.phi}, {self.phiGy}) 
      for i=1,self._grid:numCells(self._ndim) do
         self.phiGy[i]:sync()
      end
   end

   if self._isElectromagnetic then
      -- Get electromagnetic terms.
      self.apar = potentials.apar
      self.dApardt = potentials.dApardt
      self.dApardtProv = auxFields[3]
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
      self.f0Ptr = self.f0:get(1)
      self.f0lPtr = self.f0:get(1)
      self.f0Idxr = self.f0:genIndexer()

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

-- Volume integral term for use in DG scheme.
function DeltaFGyrokinetic:volTerm(w, dx, idx, f, out)
   return self.volTermFunc(w, dx, idx, f, out)
end
function DeltaFGyrokinetic:volTermGenGeo(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.f0:fill(self.f0Idxr(idx), self.f0Ptr)
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
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtProvPtr:data(), self.f0Ptr:data(), f:data(), out:data())
   else 
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.f0Ptr:data(), f:data(), out:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function DeltaFGyrokinetic:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function DeltaFGyrokinetic:surfTermGenGeo(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.f0:fill(self.f0Idxr(idxl), self.f0lPtr)
   self.f0:fill(self.f0Idxr(idxr), self.f0Ptr)
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
     res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.phiPtr:data(), self.f0lPtr:data(), self.f0Ptr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function DeltaFGyrokinetic:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f, fRefl)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.phiWall:fill(self.phiWallIdxr(idx), self.phiWallPtr)
   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, 
                                     self.phiPtr:data(), self.phiWallPtr:data(), f:data(), fRefl:data())
end

local DeltaFGyrokineticStep2 = Proto(EqBase)
-- ctor
function DeltaFGyrokineticStep2:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "DeltaFGyrokineticStep2: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "DeltaFGyrokineticStep2: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "DeltaFGyrokineticStep2: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "DeltaFGyrokineticStep2: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "DeltaFGyrokineticStep2: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   self._positivity = xsys.pickBool(tbl.positivity,false)

   self.Bvars = tbl.Bvars

   self.geoType = tbl.geometry and tbl.geometry or "SimpleHelical"

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = DeltaFGyrokineticModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p, self.geoType)
   self._surfTerm = DeltaFGyrokineticModDecl.selectStep2Surf(nm, self._cdim, self._vdim, p, self._positivity, self.Bvars, self.geoType)

   -- Select the appropriate volume and surface term functions to call.
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return DeltaFGyrokineticStep2["surfTerm" .. self.geoType](self,dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end

   self._isFirst = true
end

function DeltaFGyrokineticStep2:setAuxFields(auxFields)
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
function DeltaFGyrokineticStep2:volTerm(w, dx, idx, f, out)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.dApardt:fill(self.dApardtIdxr(idx), self.dApardtPtr)
   return self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.dApardtPtr:data(), f:data(), out:data())
end

-- Surface integral term for use in DG scheme 
-- NOTE: only vpar direction for this term
function DeltaFGyrokineticStep2:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function DeltaFGyrokineticStep2:surfTermSimpleHelical(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
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
function DeltaFGyrokineticStep2:surfTermGenGeo(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
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

return {GkEq = DeltaFGyrokinetic, GkEqStep2 = DeltaFGyrokineticStep2} 
