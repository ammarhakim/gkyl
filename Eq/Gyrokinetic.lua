-- Gkyl ------------------------------------------------------------------------
--
-- gyrokinetic equation using Hamiltonian formulation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase = require "Eq.EqBase"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local xsys = require "xsys"

local Gyrokinetic = Proto(EqBase)

-- ctor
function Gyrokinetic:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "Gyrokinetic: must specify a grid")
   self._basis = assert(tbl.phaseBasis, "Gyrokinetic: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "Gyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   -- set hamiltonian discontinuity direction flags.. CURRENTLY NOT SUPPORTED
   --self._hamilDisCont = {}
   --for d = 1, self._ndim do
   --   self._hamilDisCont[d] = false
   --end
   --if tbl.hamilDisContDirs then
   --   for i, d in ipairs(tbl.hamilDisContDirs) do
   --      self._hamilDisCont[d] = true
   --   end
   --end

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self.charge = charge
   self.mass = mass

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)

   self.Bvars = tbl.Bvars

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars)

   -- for sheath BCs
   if tbl.hasSheathBcs then
      self._calcSheathDeltaPhi = GyrokineticModDecl.selectSheathDeltaPhi(nm, self._cdim, p)
      self._calcSheathPartialReflection = GyrokineticModDecl.selectSheathPartialReflection(nm, self._cdim, self._vdim, p)
   end

   self._isFirst = true

   -- timers
   self.totalVolTime = 0.0
   self.totalSurfTime = 0.0
end

function Gyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1] -- first auxField is Field object
   local geo = auxFields[2] -- second auxField is FuncField object

   -- get phi
   self.phi = potentials.phi

   if self._isElectromagnetic then
      -- get electromagnetic terms
      self.apar = potentials.apar
      self.dApardt = potentials.dApardt
   end

   -- get magnetic geometry fields
   self.bmag = geo.bmag
   self.bmagInv = geo.bmagInv
   self.gradpar = geo.gradpar
   self.geoX = geo.geoX
   self.geoY = geo.geoY
   self.geoZ = geo.geoZ
   self.phiWall = geo.phiWall  -- for sheath BCs

   if self._isFirst then
      -- allocate pointers and indexers to field objects

      -- potentials
      self.phiPtr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      if self._isElectromagnetic then
         self.aparPtr = self.apar:get(1)
         self.dApardtPtr = self.dApardt:get(1)
         self.aparIdxr = self.apar:genIndexer()
         self.dApardtIdxr = self.dApardt:genIndexer()
      end

      -- geometry
      self.bmagPtr = self.bmag:get(1)
      self.bmagInvPtr = self.bmagInv:get(1)
      self.gradparPtr = self.gradpar:get(1)
      self.geoXPtr = self.geoX:get(1)
      self.geoYPtr = self.geoY:get(1)
      self.geoZPtr = self.geoZ:get(1)
      self.phiWallPtr = self.phiWall:get(1)
      self.bmagIdxr = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.gradparIdxr = self.gradpar:genIndexer()
      self.geoXIdxr = self.geoX:genIndexer()
      self.geoYIdxr = self.geoY:genIndexer()
      self.geoZIdxr = self.geoZ:genIndexer()
      self.phiWallIdxr = self.phiWall:genIndexer()

      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function Gyrokinetic:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idx), self.bmagInvPtr)
   self.gradpar:fill(self.gradparIdxr(idx), self.gradparPtr)
   self.geoX:fill(self.geoXIdxr(idx), self.geoXPtr)
   self.geoY:fill(self.geoYIdxr(idx), self.geoYPtr)
   self.geoZ:fill(self.geoZIdxr(idx), self.geoZPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idx), self.aparPtr)
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.geoXPtr:data(), self.geoYPtr:data(), self.geoZPtr:data(), self.phiPtr:data(), self.aparPtr:data(), f:data(), out:data())
   else 
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.geoXPtr:data(), self.geoYPtr:data(), self.geoZPtr:data(), self.phiPtr:data(), f:data(), out:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme
function Gyrokinetic:surfTerm(dir, dt, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.gradpar:fill(self.gradparIdxr(idxr), self.gradparPtr)
   self.geoX:fill(self.geoXIdxr(idxr), self.geoXPtr)
   self.geoY:fill(self.geoYIdxr(idxr), self.geoYPtr)
   self.geoZ:fill(self.geoZIdxr(idxr), self.geoZPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
     self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)
     res = self._surfTerms[dir](self.charge, self.mass, dt, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.geoXPtr:data(), self.geoYPtr:data(), self.geoZPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   else 
     res = self._surfTerms[dir](self.charge, self.mass, dt, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.geoXPtr:data(), self.geoYPtr:data(), self.geoZPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

-- calculate deltaPhi at domain edge for sheath BCs
function Gyrokinetic:calcSheathDeltaPhi(idx, edgeVal) 
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.phiWall:fill(self.phiWallIdxr(idx), self.phiWallPtr)
   return self._calcSheathDeltaPhi(self.phiPtr:data(), self.phiWallPtr:data(), edgeVal)
end
-- calculate deltaPhi at domain edge for sheath BCs
function Gyrokinetic:calcSheathPartialReflection(w, dv, edgeVal, vcut, f, fhat)
   return self._calcSheathPartialReflection(w, dv, edgeVal, vcut, f:data(), fhat:data())
end

local GyrokineticStep2 = Proto(EqBase)
-- ctor
function GyrokineticStep2:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "GyrokineticStep2: must specify a grid")
   self._basis = assert(tbl.phaseBasis, "GyrokineticStep2: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "GyrokineticStep2: must specify confBasis")

   self._ndim = self._grid:ndim()
   local charge = assert(tbl.charge, "GyrokineticStep2: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "GyrokineticStep2: must specify mass using 'mass' ")
   self.charge = charge
   self.mass = mass

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   self.Bvars = tbl.Bvars

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = GyrokineticModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, true, self.Bvars)

   self._isFirst = true
end

function GyrokineticStep2:setAuxFields(auxFields)
   local potentials = auxFields[1] -- first auxField is Field object
   local geo = auxFields[2] -- second auxField is FuncField object

   -- get phi, Apar, and dApar/dt
   self.phi = potentials.phi
   self.apar = potentials.apar
   self.dApardt = potentials.dApardt

   -- get magnetic geometry fields
   self.bmag = geo.bmag
   self.bmagInv = geo.bmagInv
   self.gradpar = geo.gradpar
   self.geoX = geo.geoX
   self.geoY = geo.geoY
   self.geoZ = geo.geoZ

   if self._isFirst then
      -- allocate pointers and indexers to field objects

      -- potentials
      self.phiPtr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      self.aparPtr = self.apar:get(1)
      self.dApardtPtr = self.dApardt:get(1)
      self.aparIdxr = self.apar:genIndexer()
      self.dApardtIdxr = self.dApardt:genIndexer()

      -- geometry
      self.bmagPtr = self.bmag:get(1)
      self.bmagInvPtr = self.bmagInv:get(1)
      self.gradparPtr = self.gradpar:get(1)
      self.geoXPtr = self.geoX:get(1)
      self.geoYPtr = self.geoY:get(1)
      self.geoZPtr = self.geoZ:get(1)
      self.bmagIdxr = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.gradparIdxr = self.gradpar:genIndexer()
      self.geoXIdxr = self.geoX:genIndexer()
      self.geoYIdxr = self.geoY:genIndexer()
      self.geoZIdxr = self.geoZ:genIndexer()

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
function GyrokineticStep2:surfTerm(dir, dt, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.gradpar:fill(self.gradparIdxr(idxr), self.gradparPtr)
   self.geoX:fill(self.geoXIdxr(idxr), self.geoXPtr)
   self.geoY:fill(self.geoYIdxr(idxr), self.geoYPtr)
   self.geoZ:fill(self.geoZIdxr(idxr), self.geoZPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)

   local res = self._surfTerms[dir](self.charge, self.mass, dt, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.geoXPtr:data(), self.geoYPtr:data(), self.geoZPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())

   return res
end

return {GkEq = Gyrokinetic, GkEqStep2 = GyrokineticStep2} 
