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
   self._nMoments = tbl.nMoments

   -- store pointers to C kernels implementing volume and surface terms
   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = PassiveAdvectionModDecl.selectVol(nm, self._ndim, p)
   self._surfTerms = PassiveAdvectionModDecl.selectSurf(nm, self._ndim, p)

   self._isFirst = true

   -- timers
   self.totalVolTime = 0.0
   self.totalSurfTime = 0.0
end

function PassiveAdvection:setAuxFields(auxFields)
end

-- Volume integral term for use in DG scheme
function PassiveAdvection:volTerm(w, dx, idx, f_ptr, fRhs_ptr)
   local tmStart = Time.clock()
   local cflRate
   cflRate = self._volTerm(w:data(), dx:data(), f_ptr:data(), fRhs_ptr:data())
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return cflRate
end

-- Surface integral term for use in DG scheme
function PassiveAdvection:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   local tmStart = Time.clock()
   local res
   res = self._surfTerms[dir](wr:data(), dxr:data(), f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data())
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

return PassiveAdvection
