-- Gkyl ------------------------------------------------------------------------
--
-- gyrokinetic equation using Hamiltonian formulation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CanonicalModDecl = require "Eq.pbData.CanonicalModDecl"
local Proto = require "Lib.Proto"
local HamiltonianBase = require "Eq.HamiltonianBase"
local ProjectOnBasis = require "Updater.ProjectOnBasis"

-- start from HamiltonianBase class
local Gyrokinetic = Proto(HamiltonianBase)

-- ctor
function Gyrokinetic:init(tbl)
   -- initialize generic Hamiltonian equation
   Gyrokinetic.super.init(self, tbl)
   assert(self._confBasis, "Gyrokinetic: must specify confBasis")

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self._qbym = charge/mass -- only q/m ratio is ever needed
   self._mcByq = mass/charge -- only q/m ratio is ever needed

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = assert(tbl.hasApar==false, "Gyrokinetic: electromagnetic not yet implemented!")

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   self._volTerm = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p)

   self._isFirst = true
end

function Gyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1]
   local funcFields = auxFields[2]

   -- get phi, and scale by q/m
   self.phi = potentials.phi
   self.phi:scale(self._qbym)

   self.bmag = funcFields.bmag
   self.bcurvY = funcFields.bcurvY

   if self._isFirst then
      -- allocate pointers to field objects
      self.phiPtr, self.phiPtrL, self.phiPtrR = self.phi:get(1), self.phi:get(1), self.phi:get(1)
      self.bmagPtr, self.bmagPtrL, self.bmagPtrR = self.bmag:get(1), self.bmag:get(1), self.bmag:get(1)
      self.bcurvYPtr, self.bcurvYPtrL, self.bcurvYPtrR = self.bcurvY:get(1), self.bcurvY:get(1), self.bcurvY:get(1)
      self.phiIdxr = self.phi:genIndexer()
      self.bmagIdxr = self.bmag:genIndexer()
      self.bcurvYIdxr = self.bcurvY:genIndexer()
      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function Gyrokinetic:volTerm(w, dx, idx, f, out)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.bcurvY:fill(self.bcurvYIdxr(idx), self.bcurvYPtr)
   return self._volTerm(self._mcByq, w:data(), dx:data(), self.bmagPtr:data(), self.bcurvYPtr:data(), self.phiPtr:data(), f:data(), out:data())
end

-- Surface integral term for use in DG scheme
function Gyrokinetic:surfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   self.phi:fill(self.phiIdxr(idxl), self.phiPtrL)
   self.phi:fill(self.phiIdxr(idxr), self.phiPtrR)
   self.bmag:fill(self.bmagIdxr(idxl), self.bmagPtrL)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtrR)
   self.bcurvY:fill(self.bcurvYIdxr(idxl), self.bcurvYPtrL)
   self.bcurvY:fill(self.bcurvYIdxr(idxr), self.bcurvYPtrR)
   return self._surfTerms[dir](self._mcByq, wr:data(), dxr:data(), maxs, self.bmagPtrR:data(), self.bcurvYPtrR:data(), self.phiPtrR:data(), fl:data(), fr:data(), outl:data(), outr:data())
end

return Gyrokinetic
