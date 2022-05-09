-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic equation using Hamiltonian formulation, except that the
-- electric field term is added outside the Poisson bracket and is in terms of
-- Epar instead of phi.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase  = require "Eq.EqBase"
local Proto   = require "Lib.Proto"
local Time    = require "Lib.Time"
local ModDecl = require "Eq.eparGkData.EparGyrokineticModDecl"
local Lin     = require "Lib.Linalg"

local EparGyrokinetic = Proto(EqBase)

-- ctor
function EparGyrokinetic:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "EparGyrokinetic: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "EparGyrokinetic: must specify a phaseBasis")
   self._confGrid  = assert(tbl.confGrid, "EparGyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "EparGyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "EparGyrokinetic: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "EparGyrokinetic: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = ModDecl.selectVol(nm, self._cdim, self._vdim, p)
   self._surfTerm = ModDecl.selectSurf(nm, self._cdim, self._vdim, p)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function EparGyrokinetic:setAuxFields(auxFields)
   local fields = auxFields[1]   -- First auxField is Field object.
   local geo    = auxFields[2]   -- Second auxField is ExternalField object.

   -- Get the parallel electric field, Epar.
   self.Epar = fields.Epar

   -- Get magnetic geometry fields.
   self.bmag        = geo.bmag
   self.cmag        = geo.cmag
   self.b_x         = geo.b_x
   self.b_y         = geo.b_y
   self.b_z         = geo.b_z
   self.jacobTotInv = geo.jacobTotInv

   if self._isFirst then
      -- Allocate pointers and indexers to field objects.

      -- Potentials.
      self.EparPtr  = self.Epar:get(1)
      self.EparIdxr = self.Epar:genIndexer()

      -- Geometry.
      self.bmagPtr         = self.bmag:get(1)
      self.cmagPtr         = self.cmag:get(1)
      self.jacobTotInvPtr  = self.jacobTotInv:get(1)
      self.b_xPtr          = self.b_x:get(1)
      self.b_yPtr          = self.b_y:get(1)
      self.b_zPtr          = self.b_z:get(1)
      self.bmagIdxr        = self.bmag:genIndexer()
      self.cmagIdxr        = self.cmag:genIndexer()
      self.jacobTotInvIdxr = self.jacobTotInv:genIndexer()
      self.b_xIdxr         = self.b_x:genIndexer()
      self.b_yIdxr         = self.b_y:genIndexer()
      self.b_zIdxr         = self.b_z:genIndexer()

      self.EparPtrFlipped  = Lin.Vec(self._confBasis:numBasis())
      self.bmagPtrFlipped  = Lin.Vec(self._confBasis:numBasis())
      self.cmagPtrFlipped  = Lin.Vec(self._confBasis:numBasis())
      self.jacobTotInvPtrFlipped  = Lin.Vec(self._confBasis:numBasis())
      self.b_xPtrFlipped          = Lin.Vec(self._confBasis:numBasis())
      self.b_yPtrFlipped          = Lin.Vec(self._confBasis:numBasis())
      self.b_zPtrFlipped          = Lin.Vec(self._confBasis:numBasis())

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function EparGyrokinetic:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.Epar:fill(self.EparIdxr(idx), self.EparPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idx), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idx), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idx), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idx), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idx), self.b_zPtr)
   local res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.EparPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function EparGyrokinetic:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   local res
   if dir==1 and idxr[1]~=1 then
      -- Use the L kernel
      self.Epar:fill(self.EparIdxr(idxl), self.EparPtr)
      self.bmag:fill(self.bmagIdxr(idxl), self.bmagPtr)
      self.jacobTotInv:fill(self.jacobTotInvIdxr(idxl), self.jacobTotInvPtr)
      self.cmag:fill(self.cmagIdxr(idxl), self.cmagPtr)
      self.b_x:fill(self.b_xIdxr(idxl), self.b_xPtr)
      self.b_y:fill(self.b_yIdxr(idxl), self.b_yPtr)
      self.b_z:fill(self.b_zIdxr(idxl), self.b_zPtr)
      res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.EparPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   else
      -- Use the R kernel
      self.Epar:fill(self.EparIdxr(idxr), self.EparPtr)
      self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
      self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
      self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
      self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
      self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
      self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
      res = self._surfTerm[dir+1](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.EparPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return {GkEq = EparGyrokinetic}
