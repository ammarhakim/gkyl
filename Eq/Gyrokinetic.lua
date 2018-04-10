-- Gkyl ------------------------------------------------------------------------
--
-- gyrokinetic equation using Hamiltonian formulation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local Proto = require "Lib.Proto"
local xsys = require "xsys"

local Gyrokinetic = Proto()

-- ctor
function Gyrokinetic:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "Gyrokinetic: must specify a grid")
   self._basis = assert(tbl.phaseBasis, "Gyrokinetic: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "Gyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   -- set hamiltonian discontinuity direction flags
   self._hamilDisCont = {}
   for d = 1, self._ndim do
      self._hamilDisCont[d] = false
   end
   if tbl.hamilDisContDirs then
      for i, d in ipairs(tbl.hamilDisContDirs) do
         self._hamilDisCont[d] = true
      end
   end

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self.charge = charge
   self.mass = mass

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)
   assert(self._isElectromagnetic==false, "Gyrokinetic: electromagnetic not yet implemented!")

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p)

   self._isFirst = true
   self.phiPtr = nil
   self.bmagPtr = nil
   self.bcurvYPtr = nil
   self.phiIdxr = nil
   self.bmagIdxr = nil
   self.bcurvYIdxr = nil
end

function Gyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1] -- first auxField is Field object
   local geo = auxFields[2] -- second auxField is FuncField object

   -- get phi
   self.phi = potentials.phi

   -- get magnetic geometry fields
   self.bmag = geo.bmag
   self.bmagInv = geo.bmagInv
   self.bcurvY = geo.bcurvY

   if self._isFirst then
      -- allocate pointers to field objects
      self.phiPtr = self.phi:get(1)
      self.bmagPtr = self.bmag:get(1)
      self.bmagInvPtr = self.bmagInv:get(1)
      self.bcurvYPtr = self.bcurvY:get(1)
      self.phiIdxr = self.phi:genIndexer()
      self.bmagIdxr = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.bcurvYIdxr = self.bcurvY:genIndexer()
      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function Gyrokinetic:volTerm(w, dx, idx, f, out)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idx), self.bmagInvPtr)
   self.bcurvY:fill(self.bcurvYIdxr(idx), self.bcurvYPtr)
   return self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.bcurvYPtr:data(), self.phiPtr:data(), f:data(), out:data())
end

-- Surface integral term for use in DG scheme
function Gyrokinetic:surfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.bcurvY:fill(self.bcurvYIdxr(idxr), self.bcurvYPtr)
   return self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.bcurvYPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
end

return Gyrokinetic
