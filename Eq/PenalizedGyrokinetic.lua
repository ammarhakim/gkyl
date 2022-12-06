-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic equation using Hamiltonian formulation with a penalization
-- factor multiplying the mu*B term in the Hamiltonian.
--
-- General geometry only kernels and no positivity are assumed.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local EqBase     = require "Eq.EqBase"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local xsys       = require "xsys"
local ModDecl    = require "Eq.penalizedGkData.PenalizedGyrokineticModDecl"
local Updater    = require "Updater"

local PenalizedGyrokinetic = Proto(EqBase)

-- ctor
function PenalizedGyrokinetic:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "PenalizedGyrokinetic: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "PenalizedGyrokinetic: must specify a phaseBasis")
   self._confGrid  = assert(tbl.confGrid, "PenalizedGyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "PenalizedGyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "PenalizedGyrokinetic: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "PenalizedGyrokinetic: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)

   self.Bvars = tbl.Bvars

   self.penalty = assert(tbl.penalty, "PenalizedGyrokinetic: must specify the penalty field in 'penalty'.")
   self.penaltyPtr  = self.penalty:get(1)
   self.penaltyIdxr = self.penalty:genIndexer()

   self._nonlinPol = xsys.pickBool(tbl.nonlinearPolarization, false)

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = ModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars)
   self._surfTerm = ModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars)

   -- Select the appropriate volume and surface term functions to call.
   local ESorEMstr = self._isElectromagnetic and "EM" or "ES"
   self.volTermFunc  = function(w, dx, idx, f, out) return PenalizedGyrokinetic["volTermGenGeo" .. ESorEMstr](self, w, dx, idx, f, out) end
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return PenalizedGyrokinetic["surfTermGenGeo" .. ESorEMstr](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
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

   self.exbEnergy = DataStruct.Field {
      onGrid        = self._confGrid,
      numComponents = self._confBasis:numBasis(),
      ghost         = {1, 1}
   }
   self.exbEnergy:clear(0.)
   if self._nonlinPol then
      -- Compute the ExB energy term in the Hamiltonian. This should go outside of the
      -- equation object but we do it here as a hack for a proof of principle first.
      -- Only available for 1x sims at the moment.
      assert(self._cdim == 1, "PenalizedGyrokinetic: nonlinear polarization only available for cdim=1.")
      self._kperpSq = assert(tbl.kperpSq, "PenalizedGyrokinetic: must specify the perpendicular wavenumber squared in 'kperpSq'.")
      self._B0 = assert(tbl.B0, "PenalizedGyrokinetic: must specify the reference magnetic field used in the polarization density with 'B0'.")
--      self._bmagInvSq = assert(tbl.bmagInvSq, "PenalizedGyrokinetic: must specify 1/B^2 with 'bmagInvSq'.")
      self.confWeakMultiply = Updater.CartFieldBinOp {
         onGrid    = self._confGrid,   operation = "Multiply",
         weakBasis = self._confBasis,  onGhosts  = true,
      }
      self.computeExBenergy = function(phiIn)
         self.confWeakMultiply:advance(0., {phiIn,phiIn}, {self.exbEnergy})
         self.exbEnergy:scale(0.5*self.mass*self._kperpSq/(self._B0^2))
--         self.confWeakMultiply:advance(0., {self._bmagInvSq,self.exbEnergy}, {self.exbEnergy})
--         self.exbEnergy:scale(0.5*self.mass*self._kperpSq)
      end
   else
      self.computeExBenergy = function(phiIn) end
   end

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function PenalizedGyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   self.computeExBenergy(self.phi) -- Compute the v_E^2 term in the Hamiltonian.

   if self._isElectromagnetic then
      -- Get electromagnetic terms.
      self.apar = potentials.apar
      self.dApardt = potentials.dApardt
      self.dApardtProv = auxFields[3]
   end

   if self._isFirst then
      -- Get magnetic geometry fields. Inside _isFirst because they are time independent.
      self.bmag = geo.bmag
      self.cmag = geo.cmag
      self.b_x = geo.b_x
      self.b_y = geo.b_y
      self.b_z = geo.b_z
      self.jacobTotInv = geo.jacobTotInv

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
      self.exbEnergyPtr  = self.exbEnergy:get(1)
      self.exbEnergyIdxr = self.exbEnergy:genIndexer()

      -- Geometry.
      self.bmagPtr     = self.bmag:get(1)
      self.cmagPtr     = self.cmag:get(1)
      self.bmagIdxr    = self.bmag:genIndexer()
      self.cmagIdxr    = self.cmag:genIndexer()
      self.jacobTotInvPtr  = self.jacobTotInv:get(1)
      self.b_xPtr         = self.b_x:get(1)
      self.b_yPtr         = self.b_y:get(1)
      self.b_zPtr         = self.b_z:get(1)
      self.jacobTotInvIdxr = self.jacobTotInv:genIndexer()
      self.b_xIdxr        = self.b_x:genIndexer()
      self.b_yIdxr        = self.b_y:genIndexer()
      self.b_zIdxr        = self.b_z:genIndexer()

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function PenalizedGyrokinetic:volTerm(w, dx, idx, f, out)
   return self.volTermFunc(w, dx, idx, f, out)
end
function PenalizedGyrokinetic:volTermGenGeoES(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.exbEnergy:fill(self.exbEnergyIdxr(idx), self.exbEnergyPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idx), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idx), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idx), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idx), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idx), self.b_zPtr)
   self.penalty:fill(self.penaltyIdxr(idx), self.penaltyPtr)
   local res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.penaltyPtr:data(), self.phiPtr:data(), self.exbEnergyPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end
function PenalizedGyrokinetic:volTermGenGeoEM(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.exbEnergy:fill(self.exbEnergyIdxr(idx), self.exbEnergyPtr)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idx), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idx), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idx), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idx), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idx), self.b_zPtr)
   self.penalty:fill(self.penaltyIdxr(idx), self.penaltyPtr)
   self.apar:fill(self.aparIdxr(idx), self.aparPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idx), self.dApardtProvPtr)
   local res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.penaltyPtr:data(), self.phiPtr:data(), self.exbEnergyPtr:data(), self.aparPtr:data(), self.dApardtProvPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function PenalizedGyrokinetic:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function PenalizedGyrokinetic:surfTermGenGeoES(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.exbEnergy:fill(self.exbEnergyIdxr(idxr), self.exbEnergyPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
   self.penalty:fill(self.penaltyIdxr(idxr), self.penaltyPtr)
   local res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.penaltyPtr:data(), self.phiPtr:data(), self.exbEnergyPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end
function PenalizedGyrokinetic:surfTermGenGeoEM(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.exbEnergy:fill(self.exbEnergyIdxr(idxr), self.exbEnergyPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
   self.penalty:fill(self.penaltyIdxr(idxr), self.penaltyPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.apar:fill(self.aparIdxr(idxl), self.aparLPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)
   self.emMod:fill(self.emModIdxr(idxl), self.emModPtrL)
   self.emMod:fill(self.emModIdxr(idxr), self.emModPtrR)
   local res = self._surfTerm[dir](self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.penaltyPtr:data(), self.phiPtr:data(), self.exbEnergyPtr:data(), self.aparPtr:data(), self.aparLPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), self.emModPtrL:data(), self.emModPtrR:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end


local PenalizedGyrokineticStep2 = Proto(EqBase)
-- ctor
function PenalizedGyrokineticStep2:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "PenalizedGyrokineticStep2: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "PenalizedGyrokineticStep2: must specify a phaseBasis")
   self._confBasis = assert(tbl.confBasis, "PenalizedGyrokineticStep2: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "PenalizedGyrokineticStep2: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "PenalizedGyrokineticStep2: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   self.penalty = assert(tbl.penalty, "PenalizedGyrokineticStep2: must specify the penalty field in 'penalty'.")
   self.penaltyPtr  = self.penalty:get(1)
   self.penaltyIdxr = self.penalty:genIndexer()

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   self.Bvars = tbl.Bvars

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = ModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p)
   self._surfTerm = ModDecl.selectStep2Surf(nm, self._cdim, self._vdim, p, self.Bvars)

   -- Select the appropriate volume and surface term functions to call.
   self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
      return PenalizedGyrokineticStep2["surfTermGenGeo"](self,dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   end

   self._isFirst = true
end

function PenalizedGyrokineticStep2:setAuxFields(auxFields)
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
   self.b_x = geo.b_x
   self.b_y = geo.b_y
   self.b_z = geo.b_z
   self.jacobTotInv = geo.jacobTotInv

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
      self.jacobTotInvPtr  = self.jacobTotInv:get(1)
      self.b_xPtr         = self.b_x:get(1)
      self.b_yPtr         = self.b_y:get(1)
      self.b_zPtr         = self.b_z:get(1)
      self.jacobTotInvIdxr = self.jacobTotInv:genIndexer()
      self.b_xIdxr        = self.b_x:genIndexer()
      self.b_yIdxr        = self.b_y:genIndexer()
      self.b_zIdxr        = self.b_z:genIndexer()

      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function PenalizedGyrokineticStep2:volTerm(w, dx, idx, f, out)
   self.bmag:fill(self.bmagIdxr(idx), self.bmagPtr)
   self.dApardt:fill(self.dApardtIdxr(idx), self.dApardtPtr)
   return self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.dApardtPtr:data(), f:data(), out:data())
end

-- Surface integral term for use in DG scheme 
-- NOTE: only vpar direction for this term
function PenalizedGyrokineticStep2:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
end
function PenalizedGyrokineticStep2:surfTermGenGeo(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.jacobTotInv:fill(self.jacobTotInvIdxr(idxr), self.jacobTotInvPtr)
   self.cmag:fill(self.cmagIdxr(idxr), self.cmagPtr)
   self.b_x:fill(self.b_xIdxr(idxr), self.b_xPtr)
   self.b_y:fill(self.b_yIdxr(idxr), self.b_yPtr)
   self.b_z:fill(self.b_zIdxr(idxr), self.b_zPtr)
   self.penalty:fill(self.penaltyIdxr(idxr), self.penaltyPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.apar:fill(self.aparIdxr(idxl), self.aparLPtr)
   self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)

   local res = self._surfTerm(self.charge, self.mass, cfll, cflr, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.jacobTotInvPtr:data(), self.cmagPtr:data(), self.b_xPtr:data(), self.b_yPtr:data(), self.b_zPtr:data(), self.penaltyPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.aparLPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), nil, nil)

   return res
end

return {GkEq = PenalizedGyrokinetic, GkEqStep2 = PenalizedGyrokineticStep2} 
