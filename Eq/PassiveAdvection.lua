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

   self.cflRateByDir = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.cflRateByDirIdxr = self.cflRateByDir:genIndexer()
end

function PassiveAdvection:setAuxFields(auxFields)
end

-- Volume integral term for use in DG scheme
function PassiveAdvection:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   local res
   local cflRateByDirPtr = self.cflRateByDir:get(1)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idx), cflRateByDirPtr)
   res = self._volTerm(w:data(), dx:data(), cflRateByDirPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme
function PassiveAdvection:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   local res
   local cflRateByDirL = self.cflRateByDir:get(1)
   local cflRateByDirR = self.cflRateByDir:get(1)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idxl), cflRateByDirL)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idxr), cflRateByDirR)
   res = self._surfTerms[dir](cflRateByDirL:data(), cflRateByDirR:data(), wr:data(), dxr:data(), dtApprox, fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function PassiveAdvection:sync()
   self.cflRateByDir:sync()
end

return PassiveAdvection
