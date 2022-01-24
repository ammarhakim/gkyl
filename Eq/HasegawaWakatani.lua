-- Gkyl ------------------------------------------------------------------------
--
-- 2D incompressible Euler equations.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local EqBase     = require "Eq.EqBase"
local ModDecl    = require "Eq.hasegawaWakataniData.hasegawa_wakatani_mod_decl"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local xsys       = require "xsys"

-- Start from HamiltonianBase class.
local HasegawaWakatani = Proto(EqBase)

-- ctor
function HasegawaWakatani:init(tbl)
   -- Fet grid and basis.
   self._grid  = assert(tbl.onGrid, "HasegawaWakatani: must specify a grid")
   self._basis = assert(tbl.basis, "HasegawaWakatani: must specify a basis")
   self._ndim  = assert(self._grid:ndim()==2, "HasegawaWakatani: grid must be 2D")
   assert(self._ndim == 2, "Incompressible Euler equations only implemented in 2D")

   self.adiabatic_C = assert(tbl.adiabaticity, "HasegawaWakatani: must specify the adiabaticity parameter, C.")
   self.gradient_kappa = assert(tbl.gradientScale, "HasegawaWakatani: must specify the normalized gradient length scale, kappa.")
   -- Store pointers to C kernels implementing volume and surface terms.
   local nm, p = self._basis:id(), self._basis:polyOrder()

   -- Use 1x,1v canonical poisson bracket.
   self._volTerm   = ModDecl.selectVol(nm, p)
   self._surfTerms = ModDecl.selectSurf(nm, p)

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0
end

function HasegawaWakatani:setAuxFields(auxFields)
   -- Get streamfunction, phi.
   self.phi = auxFields[1].phi
   if self._isFirst then
      self.phiPtr   = self.phi:get(1)
      self.phiIdxr  = self.phi:genIndexer()
      self._isFirst = false
   end
end

-- Volume integral term for use in DG scheme.
function HasegawaWakatani:volTerm(w, dx, idx, f, out)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   local res = self._volTerm(self.adiabatic_C, self.gradient_kappa, w:data(), dx:data(), self.phiPtr:data(), f:data(), out:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return res
end

-- Surface integral term for use in DG scheme.
function HasegawaWakatani:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, fl, fr, outl, outr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phiPtr)
   local res = self._surfTerms[dir](self.adiabatic_C, self.gradient_kappa, cfll, cflr, wr:data(), dxr:data(), maxs, self.phiPtr:data(), fl:data(), fr:data(), outl:data(), outr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return HasegawaWakatani
