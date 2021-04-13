-- Gkyl ------------------------------------------------------------------------
--
-- Gyrofluid model.
-- Based on the Beer+Hammett work.
-- Initially intended to complement 1D gyrokinetic mirror work. 
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
local GyrofluidModDecl = require "Eq.gyrofluidData.gyrofluid_mod_decl"
-- In order to compute dBdz, but really that should be done in the field object.
local math = require("sci.math").generic
local diff = require("sci.diff")

local Gyrofluid = Proto(EqBase)

-- ctor
function Gyrofluid:init(tbl)

   -- Get grid and basis.
   self._grid  = assert(tbl.onGrid, "Gyrofluid: must specify a grid.")
   self._basis = assert(tbl.basis, "Gyrofluid: must specify a basis.")

   self._ndim = self._grid:ndim()

   self.charge = assert(tbl.charge, "Gyrofluid: must specify charge using 'charge'.")
   self.mass   = assert(tbl.mass, "Gyrofluid: must specify mass using 'mass'.")

   self.kappaPar   = assert(tbl.kappaPar, "Gyrofluid: must specify parallel heat conductivity using 'kappaPar'.")
   self.kapparPerp = assert(tbl.kappaPerp, "Gyrofluid: must specify perpendicular heat conductivity using 'kappaPerp'.")

   self.kPerpSq  = 0.0
   self.bmagFunc = assert(tbl.bmagFunc, "Gyrofluid: must specify the function defining the magnetic field amplitude using 'bmagFunc'.")       

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm  = GyrofluidModDecl.selectVol (nm, self._ndim, p)
   self._surfTerm = GyrofluidModDecl.selectSurf(nm, self._ndim, p)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function Gyrofluid:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is ExternalField object.
   local primMom    = auxFields[3]   -- Third auxField is the primitive moments field.

   -- Get the electrostatic potential, phi.
   self.phi = potentials.phi

   -- Get magnetic geometry fields.
   self.bmag    = geo.bmag
   self.rBmag   = geo.bmagInv
   self.cmag    = geo.cmag
   self.b_x     = geo.b_x
   self.b_y     = geo.b_y
   self.b_z     = geo.b_z
   self.jacob   = geo.jacobGeo
   self.phiWall = geo.phiWall  -- For sheath BCs.

   -- Primitive moments uPar, Tpar, Tperp.
   self.primMom = primMom

   self.kperpSq = geo.kperpSq

   if self._isFirst then

      -- Compose a field with J/B. Eventually we'll remove this or put it in geo.
      local weakMult = Updater.CartFieldBinOp {
         onGrid    = self._grid,
         weakBasis = self._basis,
         operation = "Multiply",
      }
      self.jacobDbmag = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1},
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),
                          charge    = self.charge,
                          mass      = self.mass,},
      }
      weakMult:advance(0., {self.jacob,self.rBmag}, {self.jacobDbmag})

      -- Compute dBdz. Eventually we'll remove this or put it in geo.
      local dBdzFunc = function (t, xn)
         local function bmagUnpack(...)
            local xn1 = {...}
            return self.bmagFunc(0, xn1)
         end
         local deriv   = diff.derivativef(bmagUnpack, #xn)
         local xntable = {}
         for i = 1, #xn do xntable[i] = xn[i] end
         local f, dx = deriv(unpack(xntable))
         return dx
      end
      self.dBdz = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1},
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),
                          charge    = self.charge,
                          mass      = self.mass,},
      }
      local evOnNodes = Updater.EvalOnNodes {
         onGrid   = self._grid,
         basis    = self._basis,
         evaluate = dBdzFunc,
         onGhosts = true,
      }
      evOnNodes:advance(0., {}, {self.dBdz})

      -- Allocate pointers and indexers to field objects.

      self.indexer = self.bmag:genIndexer()

      -- Potentials.
      self.phiPtr  = self.phi:get(1)
      self.phiPtrl = self.phi:get(1)
      self.phiPtrr = self.phi:get(1)

      -- Geometry.
      self.bmagPtr       = self.bmag:get(1)
      self.rBmagPtr      = self.rBmag:get(1)
      self.cmagPtr       = self.cmag:get(1)
      self.b_xPtr        = self.b_x:get(1)
      self.b_yPtr        = self.b_y:get(1)
      self.b_zPtr        = self.b_z:get(1)
      self.jacobPtr      = self.jacob:get(1)
      self.phiWallPtr    = self.phiWall:get(1)
      self.jacobDbmagPtr = self.jacobDbmag:get(1)
      self.dBdzPtr       = self.dBdz:get(1)

      -- Primitive moments.
      self.primMomPtr  = self.primMom:get(1)
      self.primMomPtrl = self.primMom:get(1)
      self.primMomPtrr = self.primMom:get(1)

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function Gyrofluid:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idx), self.phiPtr)
   self.bmag:fill(self.indexer(idx), self.bmagPtr)
   self.rBmag:fill(self.indexer(idx), self.rBmagPtr)
   self.jacob:fill(self.indexer(idx), self.jacobPtr)
   self.jacobDbmag:fill(self.indexer(idx), self.jacobDbmagPtr)
   self.primMom:fill(self.indexer(idx), self.primMomPtr)
   self.dBdz:fill(self.indexer(idx), self.dBdzPtr)

   local res = self._volTerm(self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, w:data(), dx:data(), self.jacobPtr:data(), self.rBmagPtr:data(), self.jacobDbmagPtr:data(), self.dBdzPtr:data(), f:data(), self.phiPtr:data(), self.primMomPtr:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function Gyrofluid:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrl)
   self.phi:fill(self.indexer(idxr), self.phiPtrr)
   self.bmag:fill(self.indexer(idxr), self.bmagPtr)
   self.rBmag:fill(self.indexer(idxr), self.rBmagPtr)
   self.jacob:fill(self.indexer(idxr), self.jacobPtr)
   self.jacobDbmag:fill(self.indexer(idxr), self.jacobDbmagPtr)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrl)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrr)

   local res = self._surfTerm[dir](self.charge, self.mass, self.kappaPar, self.kapparPerp, self.kPerpSq, wl:data(), dxl:data(), wr:data(), dxr:data(), maxs, self.jacobPtr:data(), self.rBmagPtr:data(), self.jacobDbmagPtr:data(), fl:data(), fr:data(), self.phiPtrl:data(), self.phiPtrr:data(), self.primMomPtrl:data(), self.primMomPtrr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return Gyrofluid
