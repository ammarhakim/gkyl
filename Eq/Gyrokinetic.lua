-- Gkyl ------------------------------------------------------------------------
--
-- gyrokinetic equation using Hamiltonian formulation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local EqBase = require "Eq.EqBase"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local xsys = require "xsys"
local ffi = require "ffi"
local ffiC = ffi.C

local Gyrokinetic = Proto(EqBase)

-- ctor
function Gyrokinetic:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "Gyrokinetic: must specify a grid")
   self._basis = assert(tbl.phaseBasis, "Gyrokinetic: must specify a phaseBasis")
   self._confGrid = assert(tbl.confGrid, "Gyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "Gyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self.charge = charge
   self.mass = mass

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)
   self._positivity = xsys.pickBool(tbl.positivity,false)

   self.Bvars = tbl.Bvars

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self.Bvars)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self._positivity, self.Bvars)

   -- for sheath BCs
   if tbl.hasSheathBcs then
      self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(nm, self._cdim, self._vdim, p)
   end

   if self._isElectromagnetic then
      self.emMod = DataStruct.Field {
         onGrid = self._grid,
         numComponents = self._basis:numBasis(),
         ghost = {1, 1}
      }
      self.emModPtrL = self.emMod:get(1)
      self.emModPtrR = self.emMod:get(1)
      self.emModIdxr = self.emMod:genIndexer()
      self.emMod:clear(0.0)
   end

   -- for gyroaveraging
   self.gyavgSlvr = tbl.gyavgSlvr
   if self.gyavgSlvr then
      self._gyavg = true
      self.phiGy = {}
      local nmu = self._grid:numCells(self._ndim)
      for i = 1, nmu do
         self.phiGy[i] = DataStruct.Field {
            onGrid = self._confGrid,
            numComponents = self._confBasis:numBasis(),
            ghost = {1, 1}
         }
      end
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
   if self._gyavg then 
      self.gyavgSlvr:advance(0, {self.phi}, {self.phiGy}) 
      for i=1,self._grid:numCells(self._ndim) do
         self.phiGy[i]:sync()
      end
   end

   if self._isElectromagnetic then
      -- get electromagnetic terms
      self.apar = potentials.apar
      self.dApardt = potentials.dApardt
      self.dApardtProv = auxFields[3]
   end

   -- get magnetic geometry fields
   self.bmag = geo.bmag
   self.bmagInv = geo.bmagInv
   self.gradpar = geo.gradpar
   self.bdriftX = geo.bdriftX
   self.bdriftY = geo.bdriftY
   self.phiWall = geo.phiWall  -- for sheath BCs

   if self._isFirst then
      -- allocate pointers and indexers to field objects

      -- potentials
      self.phiPtr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      if self._isElectromagnetic then
         self.aparPtr = self.apar:get(1)
         self.dApardtPtr = self.dApardt:get(1)
         self.dApardtProvPtr = self.dApardtProv:get(1)
         self.aparIdxr = self.apar:genIndexer()
         self.dApardtIdxr = self.dApardt:genIndexer()
      end

      -- for gyroaveraging
      if self._gyavg then
         self.phiGyPtr = {}
         self.phiGyIdxr = {}
         for i=1,self._grid:numCells(self._ndim) do
            self.phiGyPtr[i] = self.phiGy[i]:get(1)
            self.phiGyIdxr[i] = self.phiGy[i]:genIndexer()
         end
      end

      -- geometry
      self.bmagPtr = self.bmag:get(1)
      self.bmagInvPtr = self.bmagInv:get(1)
      self.gradparPtr = self.gradpar:get(1)
      self.bdriftXPtr = self.bdriftX:get(1)
      self.bdriftYPtr = self.bdriftY:get(1)
      self.phiWallPtr = self.phiWall:get(1)
      self.bmagIdxr = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.gradparIdxr = self.gradpar:genIndexer()
      self.bdriftXIdxr = self.bdriftX:genIndexer()
      self.bdriftYIdxr = self.bdriftY:genIndexer()
      self.phiWallIdxr = self.phiWall:genIndexer()

      self._isFirst = false -- no longer first time
   end
end

-- Volume integral term for use in DG scheme
function Gyrokinetic:volTerm(w, dx, idx, f, out)
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
   self.gradpar:fill(self.gradparIdxr(idx), self.gradparPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idx), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idx), self.bdriftYPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idx), self.aparPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idx), self.dApardtProvPtr)
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtProvPtr:data(), f:data(), out:data())
   else 
     res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), 1., self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), f:data(), out:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme
function Gyrokinetic:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
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
   self.gradpar:fill(self.gradparIdxr(idxr), self.gradparPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftYPtr)
   local res
   if self._isElectromagnetic then
     self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
     self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)
     self.emMod:fill(self.emModIdxr(idxl), self.emModPtrL)
     self.emMod:fill(self.emModIdxr(idxr), self.emModPtrR)
     res = self._surfTerms[dir](self.charge, self.mass, cfll, cflr, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), self.emModPtrL:data(), self.emModPtrR:data())
   else 
     res = self._surfTerms[dir](self.charge, self.mass, cfll, cflr, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
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

   self._positivity = xsys.pickBool(tbl.positivity,false)
   self.Bvars = tbl.Bvars

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = GyrokineticModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p)
   self._surfTerm = GyrokineticModDecl.selectStep2Surf(nm, self._cdim, self._vdim, p, self._positivity, self.Bvars)

   self._isFirst = true
end

function GyrokineticStep2:setAuxFields(auxFields)
   local potentials = auxFields[1] -- first auxField is Field object
   local geo = auxFields[2] -- second auxField is FuncField object
   --local potentialsProv = auxFields[3]

   -- get phi, Apar, and dApar/dt
   self.phi = potentials.phi
   self.apar = potentials.apar
   self.dApardt = potentials.dApardt
   self.dApardtProv = auxFields[3]

   -- get magnetic geometry fields
   self.bmag = geo.bmag
   self.bmagInv = geo.bmagInv
   self.gradpar = geo.gradpar
   self.bdriftX = geo.bdriftX
   self.bdriftY = geo.bdriftY

   if self._isFirst then
      -- allocate pointers and indexers to field objects

      -- potentials
      self.phiPtr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      self.aparPtr = self.apar:get(1)
      self.dApardtPtr = self.dApardt:get(1)
      self.dApardtProvPtr = self.dApardtProv:get(1)
      self.aparIdxr = self.apar:genIndexer()
      self.dApardtIdxr = self.dApardt:genIndexer()

      -- geometry
      self.bmagPtr = self.bmag:get(1)
      self.bmagInvPtr = self.bmagInv:get(1)
      self.gradparPtr = self.gradpar:get(1)
      self.bdriftXPtr = self.bdriftX:get(1)
      self.bdriftYPtr = self.bdriftY:get(1)
      self.bmagIdxr = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.gradparIdxr = self.gradpar:genIndexer()
      self.bdriftXIdxr = self.bdriftX:genIndexer()
      self.bdriftYIdxr = self.bdriftY:genIndexer()

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
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.bmag:fill(self.bmagIdxr(idxr), self.bmagPtr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInvPtr)
   self.gradpar:fill(self.gradparIdxr(idxr), self.gradparPtr)
   self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftXPtr)
   self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftYPtr)
   self.apar:fill(self.aparIdxr(idxr), self.aparPtr)
   self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardtPtr)
   self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProvPtr)

   local res = self._surfTerm(self.charge, self.mass, cfll, cflr, wr:data(), dxr:data(), maxs, self.bmagPtr:data(), self.bmagInvPtr:data(), self.gradparPtr:data(), self.bdriftXPtr:data(), self.bdriftYPtr:data(), self.phiPtr:data(), self.aparPtr:data(), self.dApardtPtr:data(), self.dApardtProvPtr:data(), fl:data(), fr:data(), outl:data(), outr:data(), nil, nil)

   return res
end

return {GkEq = Gyrokinetic, GkEqStep2 = GyrokineticStep2} 
