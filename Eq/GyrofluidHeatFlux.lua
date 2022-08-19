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
   print("setting kPerpSq = 0.0 in Eq.GyrofluidHeatFlux.lua.")

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = GFheatFluxModDecl.selectVol (nm, self._ndim, p)
   self._surfTerm = GFheatFluxModDecl.selectSurf(nm, self._ndim, p)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function GyrofluidHeatFlux:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.
   local primMom    = auxFields[3]   -- Third auxField is the primitive moments field.

   -- In case there is no Jacobian:
   if self._isFirst then
      if geo.jacobGeo == nil then
         self.unitFld = DataStruct.Field {
            onGrid        = self._grid,
            numComponents = self._basis:numBasis(),
            ghost         = {1, 1},
            metaData      = {polyOrder = self._basis:polyOrder(),
                             basisType = self._basis:id(),},
         }
         local evOnNodes = Updater.EvalOnNodes {
            onGrid = self._grid,   evaluate = function(t, xn) return 1. end,
            basis  = self._basis,  onGhosts = true,
         }
         evOnNodes:advance(0., {}, {self.unitFld})
      end
   end

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   -- Get magnetic geometry fields.
   self.rBmag   = geo.bmagInv
   self.jacob   = geo.jacobGeo or self.unitFld
   self.phiWall = geo.phiWall  -- For sheath BCs.

   -- Primitive moments uPar, Tpar, Tperp.
   self.primMom = primMom

   self.kperpSq = geo.kperpSq

   if self._isFirst then

      -- Compose a field with 1/B^2. Eventually we'll remove this or put it in geo.
      local weakMult = Updater.CartFieldBinOp {
         weakBasis = self._basis,  operation = "Multiply",
         onGhosts  = true,
      }
      self.rBmagSq = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1},
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),},
      }
      self.rBmagSq:clear(0.)
      weakMult:advance(0., {self.rBmag,self.rBmag}, {self.rBmagSq})

      -- Allocate pointers and indexers to field objects.

      self.indexer = self.rBmag:genIndexer()

      -- Potentials.
      self.phiPtr = self.phi:get(1)
      self.phiPtrL, self.phiPtrR = self.phi:get(1), self.phi:get(1)

      -- Geometry.
      self.jacobPtr = self.jacob:get(1)
      self.jacobPtrL, self.jacobPtrR = self.jacob:get(1), self.jacob:get(1)
      self.rBmagPtr = self.rBmag:get(1)
      self.rBmagPtrL, self.rBmagPtrR = self.rBmag:get(1), self.rBmag:get(1)
      self.rBmagSqPtr = self.rBmagSq:get(1)
      self.rBmagSqPtrL, self.rBmagSqPtrR = self.rBmagSq:get(1), self.rBmagSq:get(1)
      self.phiWallPtr = self.phiWall:get(1)

      -- Primitive moments.
      self.primMomPtr = self.primMom:get(1)
      self.primMomPtrL, self.primMomPtrR = self.primMom:get(1), self.primMom:get(1)

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function GyrofluidHeatFlux:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idx), self.phiPtr)
   self.jacob:fill(self.indexer(idx), self.jacobPtr)
   self.rBmag:fill(self.indexer(idx), self.rBmagPtr)
   self.rBmagSq:fill(self.indexer(idx), self.rBmagSqPtr)
   self.primMom:fill(self.indexer(idx), self.primMomPtr)

   local res = self._volTerm(self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, w:data(), dx:data(), self.jacobPtr:data(), self.rBmagPtr:data(), self.rBmagSqPtr:data(), f:data(), self.phiPtr:data(), self.primMomPtr:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function GyrofluidHeatFlux:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrL)
   self.rBmag:fill(self.indexer(idxl), self.rBmagPtrL)
   self.rBmagSq:fill(self.indexer(idxl), self.rBmagSqPtrL)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrL)

   self.phi:fill(self.indexer(idxr), self.phiPtrR)
   self.rBmag:fill(self.indexer(idxr), self.rBmagPtrR)
   self.rBmagSq:fill(self.indexer(idxr), self.rBmagSqPtrR)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrR)

   self._surfTerm[dir](self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, wl:data(), wr:data(), dxl:data(), dxr:data(), self.rBmagPtrL:data(), self.rBmagPtrR:data(), self.rBmagSqPtrL:data(), self.rBmagSqPtrR:data(), fl:data(), fr:data(), self.phiPtrL:data(), self.phiPtrR:data(), self.primMomPtrL:data(), self.primMomPtrR:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return 0
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function GyrofluidHeatFlux:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrL)
   self.rBmag:fill(self.indexer(idxl), self.rBmagPtrL)
   self.rBmagSq:fill(self.indexer(idxl), self.rBmagSqPtrL)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrL)

   self.phi:fill(self.indexer(idxr), self.phiPtrR)
   self.rBmag:fill(self.indexer(idxr), self.rBmagPtrR)
   self.rBmagSq:fill(self.indexer(idxr), self.rBmagSqPtrR)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrR)

   self._boundarySurfTerms[dir](self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, wl:data(), wr:data(), dxl:data(), dxr:data(), maxs, self.rBmagPtrL:data(), self.rBmagPtrR:data(), self.rBmagSqPtrL:data(), self.rBmagSqPtrR:data(), fl:data(), fr:data(), self.phiPtrL:data(), self.phiPtrR:data(), self.primMomPtrL:data(), self.primMomPtrR:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return 0
end

return GyrofluidHeatFlux
