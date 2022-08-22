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
   local cRusanov   = auxFields[4]   -- Fourth auxField is the speed used in the Rusanov numerical flux.

   if self._isFirst then
      -- In case there is no Jacobian:
      if geo.jacobGeoInv == nil then
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
   self.rJacob  = geo.jacobGeoInv or self.unitFld
   self.phiWall = geo.phiWall  -- For sheath BCs.

   self.primMom = primMom   -- Primitive moments uPar, Tpar, Tperp.
   self.cRusanov = cRusanov

   if self._isFirst then

      -- Compose a field with J/B. Eventually we'll remove this or put it in geo.
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
                          basisType = self._basis:id(),},
      }
      self.dBdz:clear(0.)
      local evOnNodes = Updater.EvalOnNodes {
         onGrid = self._grid,   evaluate = dBdzFunc,
         basis  = self._basis,  onGhosts = true,
      }
      evOnNodes:advance(0., {}, {self.dBdz})

      -- Allocate pointers and indexers to field objects.

      self.indexer = self.rBmag:genIndexer()

      -- Potentials.
      self.phiPtr = self.phi:get(1)
      self.phiPtrL, self.phiPtrR = self.phi:get(1), self.phi:get(1)

      -- Geometry.
      self.rJacobPtr  = self.rJacob:get(1)
      self.rJacobPtrL, self.rJacobPtrR = self.rJacob:get(1), self.rJacob:get(1)
      self.rBmagPtr   = self.rBmag:get(1)
      self.rBmagPtrL, self.rBmagPtrR = self.rBmag:get(1), self.rBmag:get(1)
      self.rBmagSqPtr = self.rBmagSq:get(1)
      self.rBmagSqPtrL, self.rBmagSqPtrR = self.rBmagSq:get(1), self.rBmagSq:get(1)
      self.phiWallPtr = self.phiWall:get(1)
      self.dBdzPtr    = self.dBdz:get(1)

      -- Primitive moments.
      self.primMomPtr = self.primMom:get(1)
      self.primMomPtrL, self.primMomPtrR = self.primMom:get(1), self.primMom:get(1)

      -- Speed in Rusanov numerical flux.
      self.cRusanovPtr = self.cRusanov:get(1)
      self.cRusanovPtrL, self.cRusanovPtrR = self.cRusanov:get(1), self.cRusanov:get(1)

      self._isFirst = false -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function Gyrofluid:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idx), self.phiPtr)
   self.rJacob:fill(self.indexer(idx), self.rJacobPtr)
   self.rBmag:fill(self.indexer(idx), self.rBmagPtr)
   self.primMom:fill(self.indexer(idx), self.primMomPtr)
   self.dBdz:fill(self.indexer(idx), self.dBdzPtr)
   self.cRusanov:fill(self.indexer(idx), self.cRusanovPtr)

   local res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.rJacobPtr:data(), self.rBmagPtr:data(), self.dBdzPtr:data(), f:data(), self.phiPtr:data(), self.primMomPtr:data(), self.cRusanovPtr:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function Gyrofluid:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()

   self.phi:fill(self.indexer(idxl), self.phiPtrL)
   self.rJacob:fill(self.indexer(idxl), self.rJacobPtrL)
   self.rBmag:fill(self.indexer(idxl), self.rBmagPtrL)
   self.rBmagSq:fill(self.indexer(idxl), self.rBmagSqPtrL)
   self.primMom:fill(self.indexer(idxl), self.primMomPtrL)
   self.cRusanov:fill(self.indexer(idxl), self.cRusanovPtrL)

   self.phi:fill(self.indexer(idxr), self.phiPtrR)
   self.rJacob:fill(self.indexer(idxr), self.rJacobPtrR)
   self.rBmag:fill(self.indexer(idxr), self.rBmagPtrR)
   self.rBmagSq:fill(self.indexer(idxr), self.rBmagSqPtrR)
   self.primMom:fill(self.indexer(idxr), self.primMomPtrR)
   self.cRusanov:fill(self.indexer(idxr), self.cRusanovPtrR)
   
   local res = self._surfTerm[dir](self.charge, self.mass, wl:data(), wr:data(), dxl:data(), dxr:data(), maxs, self.rJacobPtrL:data(), self.rJacobPtrR:data(), self.rBmagPtrL:data(), self.rBmagPtrR:data(), self.rBmagSqPtrL:data(), self.rBmagSqPtrR:data(), fl:data(), fr:data(), self.phiPtrL:data(), self.phiPtrR:data(), self.primMomPtrL:data(), self.primMomPtrR:data(), self.cRusanovPtrL:data(), self.cRusanovPtrR:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return Gyrofluid
