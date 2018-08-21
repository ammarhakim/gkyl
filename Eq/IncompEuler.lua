-- Gkyl ------------------------------------------------------------------------
--
-- 2D incompressible Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase = require "Eq.EqBase"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local xsys = require "xsys"

-- start from HamiltonianBase class
local IncompEuler = Proto(EqBase)

-- ctor
function IncompEuler:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "IncompEuler: must specify a grid")
   self._basis = assert(tbl.basis, "IncompEuler: must specify a basis")
   self._ndim = self._grid:ndim()
   assert(self._ndim == 2, "Incompressible Euler equations only implemented in 2D")

   self.charge = tbl.charge and tbl.charge or 1
   self.mass = tbl.mass and tbl.mass or 1

   self._positivity = xsys.pickBool(tbl.positivity,false)

   -- store pointers to C kernels implementing volume and surface terms
   local nm, p = self._basis:id(), self._basis:polyOrder()
   -- use 1x,1v canonical poisson bracket
   self._volTerm = GyrokineticModDecl.selectVol(nm, 2, 0, p, false, {0})
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, 2, 0, p, false, self._positivity, {0})

   -- set up some dummy fields so that we can use GK kernels
   self.unitField = DataStruct.Field {
        onGrid = self._grid,
        numComponents = self._basis:numBasis(),
        ghost = {1, 1},
   }
   local initUnit = Updater.ProjectOnBasis {
      onGrid = self._grid,
      basis = self._basis,
      evaluate = function (t,xn)
                    return -1.0 -- sign convention
                 end,
      projectOnGhosts = true,
   }
   initUnit:advance(0.,0.,{},{self.unitField})

   self.zeroField = DataStruct.Field {
        onGrid = self._grid,
        numComponents = self._basis:numBasis(),
        ghost = {1, 1},
   }
   self.zeroField:clear(0.0)

   self.unitPtr = self.unitField:get(1)
   self.unitIdxr = self.unitField:genIndexer()
   self.zeroPtr = self.zeroField:get(1)
   self.zeroIdxr = self.zeroField:genIndexer()

   self._isFirst = true

   -- timers
   self.totalVolTime = 0.0
   self.totalSurfTime = 0.0
end

function IncompEuler:setAuxFields(auxFields)
   -- get streamfunction, phi
   self.phi = auxFields[1].phi
   if self._isFirst then
      self.phiPtr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      self._isFirst = false
   end
end

-- Volume integral term for use in DG scheme
function IncompEuler:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.unitField:fill(self.unitIdxr(idx), self.unitPtr)
   self.zeroField:fill(self.zeroIdxr(idx), self.zeroPtr)
   local res
   res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.unitPtr:data(), self.unitPtr:data(), self.zeroPtr:data(), self.zeroPtr:data(), self.zeroPtr:data(), self.phiPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme
function IncompEuler:surfTerm(dir, cfl, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   self.unitField:fill(self.unitIdxr(idxr), self.unitPtr)
   self.zeroField:fill(self.zeroIdxr(idxr), self.zeroPtr)
   local res
   res = self._surfTerms[dir](self.charge, self.mass, cfl, wr:data(), dxr:data(), maxs, self.unitPtr:data(), self.unitPtr:data(), self.zeroPtr:data(), self.zeroPtr:data(), self.zeroPtr:data(), self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return IncompEuler
