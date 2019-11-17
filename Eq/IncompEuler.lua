-- Gkyl ------------------------------------------------------------------------
--
-- 2D incompressible Euler equations.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase             = require "Eq.EqBase"
local IncompEulerModDecl = require "Eq.incompEulerData.IncompEulerModDecl"
local DataStruct         = require "DataStruct"
local Updater            = require "Updater"
local Proto              = require "Lib.Proto"
local Time               = require "Lib.Time"
local xsys               = require "xsys"

-- Start from HamiltonianBase class.
local IncompEuler = Proto(EqBase)

-- ctor
function IncompEuler:init(tbl)
   -- Fet grid and basis.
   self._grid  = assert(tbl.onGrid, "IncompEuler: must specify a grid")
   self._basis = assert(tbl.basis, "IncompEuler: must specify a basis")
   self._ndim  = self._grid:ndim()
   assert(self._ndim == 2, "Incompressible Euler equations only implemented in 2D")

   self.charge = tbl.charge and tbl.charge or 1
   self.mass   = tbl.mass and tbl.mass or 1

   self._positivity = xsys.pickBool(tbl.positivity,false)

   -- Store pointers to C kernels implementing volume and surface terms.
   local nm, p = self._basis:id(), self._basis:polyOrder()
   -- Use 1x,1v canonical poisson bracket.
   self._volTerm   = IncompEulerModDecl.selectVol(nm, 2, p)
   self._surfTerms = IncompEulerModDecl.selectSurf(nm, 2, p, self._positivity)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0

   self.cflRateByDir = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.cflRateByDirIdxr = self.cflRateByDir:genIndexer()
end

function IncompEuler:setAuxFields(auxFields)
   -- Get streamfunction, phi.
   self.phi = auxFields[1].phi
   if self._isFirst then
      self.phiPtr   = self.phi:get(1)
      self.phiIdxr  = self.phi:genIndexer()
      self._isFirst = false
   end
end

-- Volume integral term for use in DG scheme.
function IncompEuler:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)

   local cflRateByDirPtr = self.cflRateByDir:get(1)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idx), cflRateByDirPtr)

   local res
   res = self._volTerm(self.charge, self.mass, w:data(), dx:data(), cflRateByDirPtr:data(), self.phiPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function IncompEuler:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)

   local cflRateByDirL = self.cflRateByDir:get(1)
   local cflRateByDirR = self.cflRateByDir:get(1)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idxl), cflRateByDirL)
   self.cflRateByDir:fill(self.cflRateByDirIdxr(idxr), cflRateByDirR)

   local res
   res = self._surfTerms[dir](self.charge, self.mass, cflRateByDirL:data(), cflRateByDirR:data(), wr:data(), dxr:data(), dtApprox, self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function IncompEuler:sync()
   self.cflRateByDir:sync()
end

return IncompEuler
