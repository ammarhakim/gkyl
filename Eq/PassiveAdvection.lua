-- Gkyl ------------------------------------------------------------------------
--
-- 2D incompressible Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase = require "Eq.EqBase"
local PassiveAdvectionModDecl = require "Eq.passiveAdvectionData.PassiveAdvectionModDecl"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local xsys = require "xsys"

-- start from HamiltonianBase class
local PassiveAdvection = Proto(EqBase)

-- ctor
function PassiveAdvection:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "PassiveAdvection: must specify a grid")
   self._basis = assert(tbl.basis, "PassiveAdvection: must specify a basis")
   self._ndim = self._grid:ndim()

   self._positivity = xsys.pickBool(tbl.positivity, false)

   -- store pointers to C kernels implementing volume and surface terms
   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = PassiveAdvectionModDecl.selectVol(nm, self._ndim, p)
   self._surfTerms = PassiveAdvectionModDecl.selectSurf(nm, self._ndim, p, self._positivity)

   self._isFirst = true

   -- timers
   self.totalVolTime = 0.0
   self.totalSurfTime = 0.0

   self.cflRateCtrl = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.cflRateCtrlIdxr = self.cflRateCtrl:genIndexer()
end

function PassiveAdvection:setAuxFields(auxFields)
end

-- Volume integral term for use in DG scheme
function PassiveAdvection:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   local res
   local cflRateCtrlPtr = self.cflRateCtrl:get(1)
   self.cflRateCtrl:fill(self.cflRateCtrlIdxr(idx), cflRateCtrlPtr)
   res = self._volTerm(w:data(), dx:data(), cflRateCtrlPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme
function PassiveAdvection:surfTerm(dir, cflL, cflR, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   local res
   local cflRateCtrlL = self.cflRateCtrl:get(1)
   local cflRateCtrlR = self.cflRateCtrl:get(1)
   self.cflRateCtrl:fill(self.cflRateCtrlIdxr(idxl), cflRateCtrlL)
   self.cflRateCtrl:fill(self.cflRateCtrlIdxr(idxr), cflRateCtrlR)
   res = self._surfTerms[dir](cflRateCtrlL:data(), cflRateCtrlR:data(), wr:data(), dxr:data(), maxs, fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function PassiveAdvection:sync()
   self.cflRateCtrl:sync()
end

return PassiveAdvection
