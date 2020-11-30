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

   if self._positivity then
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self._grid,
         basis  = self._basis,
      }

      self.fRhsVol = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1}
      }
      self.fRhsVol_ptr = self.fRhsVol:get(1)

      self.fRhsSurf = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
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

function IncompEuler:setAuxFields(auxFields)
   -- Get streamfunction, phi.
   self.phi = auxFields[1].phi
   if self._isFirst then
      self.phi_ptr   = self.phi:get(1)
      self.phiIdxr  = self.phi:genIndexer()
      self._isFirst = false
   end
end

function IncompEuler:volTerm(w, dx, idx, f_ptr, fRhs_ptr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idx), self.phi_ptr)
   local cflRate
   if self._positivity then
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)
      self.fRhsVol:fill(self.fRhsIdxr(idx), self.fRhsVol_ptr)
      cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self.phi_ptr:data(), f_ptr:data(), self.fRhsVol_ptr:data())
   else
      cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), self.positivityWeightByDir_ptr:data(), self.phi_ptr:data(), f_ptr:data(), fRhs_ptr:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return cflRate
end

-- Surface integral term for use in DG scheme
function IncompEuler:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   local tmStart = Time.clock()
   self.phi:fill(self.phiIdxr(idxr), self.phi_ptr)
   local res
   if self._positivity then
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxl), self.positivityWeightByDir_L_ptr)
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxr), self.positivityWeightByDir_R_ptr)
      self.fRhsSurf:fill(self.fRhsIdxr(idxl), self.fRhsSurf_L_ptr)
      self.fRhsSurf:fill(self.fRhsIdxr(idxr), self.fRhsSurf_R_ptr)

      res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), dtApprox, self.positivityWeightByDir_L_ptr:data(), self.positivityWeightByDir_R_ptr:data(),
                                 self.phi_ptr:data(), f_L_ptr:data(), f_R_ptr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
   else
      res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), self.phi_ptr:data(), f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data())
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function IncompEuler:clearRhsTerms()
   if self._positivity then
      self.fRhsVol:clear(0.0)
      self.fRhsSurf:clear(0.0)
   end
end

-- When using positivity algorithm, different parts of RHS are stored separately.
-- here we combine the parts, with some rescaling of the volume term.
function IncompEuler:getPositivityRhs(tCurr, dtApprox, fIn, fRhs)
   weightDirs = {}
   for d = 1, self._ndim do
      weightDirs[d] = d
   end
   -- fIn + fac*dt*fVol + dt*fSurf > 0.
   -- Rescale volume term by fac, and add to surface term fRhs = fRhsSurf.
   self.posRescaler:rescaleVolTerm(tCurr, dtApprox, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurf, self.fRhsVol)

   fRhs:accumulate(1.0, self.fRhsSurf, 1.0, self.fRhsVol)
end

function IncompEuler:setPositivityWeights(cflRateByCell)
   -- set total weight = positivityWeightByDir[0] = cflRateByCell in each cell
   -- other elements will be overwritten in kernels
   self.positivityWeightByDir:clear(1.0)
   self.positivityWeightByDir:scaleByCell(cflRateByCell)
end

function IncompEuler:sync()
   self.positivityWeightByDir:sync()
end

return IncompEuler
