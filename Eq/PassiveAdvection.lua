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

   self._positivity = xsys.pickBool(tbl.positivity, false)

   -- store pointers to C kernels implementing volume and surface terms
   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = PassiveAdvectionModDecl.selectVol(nm, self._ndim, p)
   self._surfTerms = PassiveAdvectionModDecl.selectSurf(nm, self._ndim, p, self._positivity)

   self._isFirst = true

   -- timers
   self.totalVolTime = 0.0
   self.totalSurfTime = 0.0

   if self._positivity then
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self._grid,
         basis  = self._basis,
      }

      self.fRhsVol = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis()*self._nMoments,
         ghost         = {1, 1}
      }
      self.fRhsVol_ptr = self.fRhsVol:get(1)

      self.fRhsSurf = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis()*self._nMoments,
         ghost         = {1, 1}
      }
      self.fRhsSurf_ptr = self.fRhsSurf:get(1)
      self.fRhsSurf_L_ptr = self.fRhsSurf:get(1)
      self.fRhsSurf_R_ptr = self.fRhsSurf:get(1)

      self.fRhsIdxr = self.fRhsVol:genIndexer()
   end

   self.positivityWeightByDir = DataStruct.Field {
      onGrid        = self._grid,
      numComponents = self._basis:numBasis(),
      ghost         = {1, 1},
   }
   self.positivityWeightByDirIdxr = self.positivityWeightByDir:genIndexer()
   self.positivityWeightByDir_ptr = self.positivityWeightByDir:get(1)
   self.positivityWeightByDir_L_ptr = self.positivityWeightByDir:get(1)
   self.positivityWeightByDir_R_ptr = self.positivityWeightByDir:get(1)
end

function PassiveAdvection:setAuxFields(auxFields)
end

-- Volume integral term for use in DG scheme
function PassiveAdvection:volTerm(w, dx, idx, f_ptr, fRhs_ptr)
   local tmStart = Time.clock()
   local cflRate
   if self._positivity then
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)
      self.fRhsVol:fill(self.fRhsIdxr(idx), self.fRhsVol_ptr)
      cflRate = self._volTerm(w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
   else
      cflRate = self._volTerm(w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), f_ptr:data(), fRhs_ptr:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return cflRate
end

-- Surface integral term for use in DG scheme
function PassiveAdvection:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   local tmStart = Time.clock()
   local res
   if self._positivity then
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxl), self.positivityWeightByDir_L_ptr)
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxr), self.positivityWeightByDir_R_ptr)
      self.fRhsSurf:fill(self.fRhsIdxr(idxl), self.fRhsSurf_L_ptr)
      self.fRhsSurf:fill(self.fRhsIdxr(idxr), self.fRhsSurf_R_ptr)

      res = self._surfTerms[dir](wr:data(), dxr:data(), dtApprox, self.positivityWeightByDir_L_ptr:data(), self.positivityWeightByDir_R_ptr:data(),
                                 f_L_ptr:data(), f_R_ptr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
   else
      res = self._surfTerms[dir](wr:data(), dxr:data(), f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function PassiveAdvection:clearRhsTerms()
   if self._positivity then
      self.fRhsVol:clear(0.0)
      self.fRhsSurf:clear(0.0)
   end
end

-- When using positivity algorithm, different parts of RHS are stored separately.
-- here we combine the parts, with some rescaling of the volume term.
function PassiveAdvection:getPositivityRhs(tCurr, dtApprox, fIn, fRhs)
   weightDirs = {}
   for d = 1, self._ndim do
      weightDirs[d] = d
   end
   -- fIn + fac*dt*fVol + dt*fSurf > 0.
   -- Rescale volume term by fac, and add to surface term fRhs = fRhsSurf.
   if dtApprox > 0 then self.posRescaler:rescaleVolTerm(tCurr, dtApprox*1.05, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurf, self.fRhsVol) end

   fRhs:accumulate(1.0, self.fRhsSurf, 1.0, self.fRhsVol)
end

function PassiveAdvection:setPositivityWeights(cflRateByCell)
   -- set total weight = positivityWeightByDir[0] = cflRateByCell in each cell
   -- other elements will be overwritten in kernels
   self.positivityWeightByDir:clear(1.0)
   self.positivityWeightByDir:scaleByCell(cflRateByCell)
end

function PassiveAdvection:sync()
   self.positivityWeightByDir:sync()
end

return PassiveAdvection
