-- Gkyl ------------------------------------------------------------------------
--
-- Gyrofluid heatflux terms based on the Beer+Hammett work.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Updater    = require "Updater"
local DataStruct = require "DataStruct"
local EqBase     = require "Eq.EqBase"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local xsys       = require "xsys"
local GFheatFluxModDecl = require "Eq.gfHeatFluxData.gyrofluid_heatflux_mod_decl"

local GyrofluidHeatFlux = Proto(EqBase)

-- ctor
function GyrofluidHeatFlux:init(tbl)

   -- Get grid and basis.
   self._grid  = assert(tbl.onGrid, "Eq.GyrofluidHeatFlux: must specify a grid.")
   self._basis = assert(tbl.basis, "Eq.GyrofluidHeatFlux: must specify a basis.")

   self._ndim = self._grid:ndim()

   self.charge = assert(tbl.charge, "Eq.GyrofluidHeatFlux: must specify charge using 'charge'.")
   self.mass   = assert(tbl.mass, "Eq.GyrofluidHeatFlux: must specify mass using 'mass'.")

   self.kappaPar   = assert(tbl.kappaPar, "Eq.GyrofluidHeatFlux: must specify parallel heat conductivity using 'kappaPar'.")
   self.kapparPerp = assert(tbl.kappaPerp, "Eq.GyrofluidHeatFlux: must specify perpendicular heat conductivity using 'kappaPerp'.")

   self.kPerpSq  = 0.
   print("setting kPerpSq = 0.0 in Eq/Gyrofluid.lua.")

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = GFheatFluxModDecl.selectVol (nm, self._ndim, p)
   self._surfTerm = GFheatFluxModDecl.selectSurf(nm, self._ndim, p)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function GyrofluidHeatFlux:setAuxFields(auxFields) end
function GyrofluidHeatFlux:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.
   local primMom    = auxFields[3]   -- Third auxField is the primitive moments field.

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   -- Get magnetic geometry fields.
   self.bmag    = geo.bmag
   self.jacob   = geo.jacobGeo
   self.phiWall = geo.phiWall  -- For sheath BCs.

   -- Primitive moments uPar, Tpar, Tperp.
   self.primMom = primMom

   self.kperpSq = geo.kperpSq

   if self._isFirst then

      -- Compose a field with J/B. Eventually we'll remove this or put it in geo.
      local weakMult = Updater.CartFieldBinOp {
         onGrid    = self._grid,   operation = "Multiply",
         weakBasis = self._basis,  onGhosts  = true,
      }
      self.jacobDbmag = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1},
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),},
      }
      self.jacobDbmag:clear(0.)
      local rBmag = geo.bmagInv
      weakMult:advance(0., {self.jacob,rBmag}, {self.jacobDbmag})

      -- Allocate pointers and indexers to field objects.

      self.indexer = self.bmag:genIndexer()

      -- Potentials.
      self.phiPtr  = self.phi:get(1)
      self.phiPtrl = self.phi:get(1)
      self.phiPtrr = self.phi:get(1)

      -- Geometry.
      self.bmagPtr       = self.bmag:get(1)
      self.jacobPtr      = self.jacob:get(1)
      self.phiWallPtr    = self.phiWall:get(1)
      self.jacobDbmagPtr = self.jacobDbmag:get(1)

      -- Primitive moments.
      self.primMomPtr  = self.primMom:get(1)
      self.primMomPtrl = self.primMom:get(1)
      self.primMomPtrr = self.primMom:get(1)

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function GyrofluidHeatFlux:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idx), self.phiPtr)
   self.jacob:fill(self.indexer(idx), self.jacobPtr)
   self.jacobDbmag:fill(self.indexer(idx), self.jacobDbmagPtr)
   self.primMom:fill(self.indexer(idx), self.primMomPtr)

   local res = self._volTerm(self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, w:data(), dx:data(), self.jacobPtr:data(), self.jacobDbmagPtr:data(), f:data(), self.phiPtr:data(), self.primMomPtr:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function GyrofluidHeatFlux:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrl)
   self.phi:fill(self.indexer(idxr), self.phiPtrr)
   self.jacob:fill(self.indexer(idxl), self.jacobPtr)
   self.jacobDbmag:fill(self.indexer(idxl), self.jacobDbmagPtr)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrl)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrr)

   local res = self._surfTerm[dir](self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, wl:data(), dxl:data(), wr:data(), dxr:data(), self.jacobPtr:data(), self.jacobDbmagPtr:data(), fl:data(), fr:data(), self.phiPtrl:data(), self.phiPtrr:data(), self.primMomPtrl:data(), self.primMomPtrr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function GyrofluidHeatFlux:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrl)
   self.phi:fill(self.indexer(idxr), self.phiPtrr)
   self.jacob:fill(self.indexer(idxl), self.jacobPtr)
   self.jacobDbmag:fill(self.indexer(idxl), self.jacobDbmagPtr)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrl)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrr)

   self._boundarySurfTerms[dir](self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.jacobPtr:data(), self.jacobDbmagPtr:data(), fl:data(), fr:data(), self.phiPtrl:data(), self.phiPtrr:data(), self.primMomPtrl:data(), self.primMomPtrr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return GyrofluidHeatFlux
