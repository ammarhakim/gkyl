-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic equation using Hamiltonian formulation.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct         = require "DataStruct"
local EqBase             = require "Eq.EqBase"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"
local Proto              = require "Lib.Proto"
local Time               = require "Lib.Time"
local Updater            = require "Updater"
local xsys               = require "xsys"
local ffi                = require "ffi"
local ffiC               = ffi.C

local Gyrokinetic = Proto(EqBase)

-- ctor
function Gyrokinetic:init(tbl)
   -- Get grid and basis.
   self._grid      = assert(tbl.onGrid, "Gyrokinetic: must specify a grid")
   self._basis     = assert(tbl.phaseBasis, "Gyrokinetic: must specify a phaseBasis")
   self._confGrid  = assert(tbl.confGrid, "Gyrokinetic: must specify confGrid")
   self._confBasis = assert(tbl.confBasis, "Gyrokinetic: must specify confBasis")

   self._ndim = self._grid:ndim()

   local charge = assert(tbl.charge, "Gyrokinetic: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "Gyrokinetic: must specify mass using 'mass' ")
   self.charge  = charge
   self.mass    = mass

   assert(tbl.hasPhi==true, "Gyrokinetic: must have an electrostatic potential!")
   self._isElectromagnetic = xsys.pickBool(tbl.hasApar, false)
   self._positivity        = xsys.pickBool(tbl.positivity,false)

   self.Bvars = tbl.Bvars

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm   = GyrokineticModDecl.selectVol(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self._positivity, self.Bvars)
   self._surfTerms = GyrokineticModDecl.selectSurf(nm, self._cdim, self._vdim, p, self._isElectromagnetic, self._positivity, self.Bvars)
   if self._isElectromagnetic then 
      self._volTermStep2 = GyrokineticModDecl.selectStep2Vol(nm, self._cdim, self._vdim, p, self._positivity)
      if p == 1 and not self._positivity then self._surfTermsStep2 = GyrokineticModDecl.selectStep2Surf(nm, self._cdim, self._vdim, p, self._positivity, self.Bvars) end
   end
   self._polyOrder = p

   -- For sheath BCs.
   if tbl.hasSheathBcs then
      self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(nm, self._cdim, self._vdim, p)
   end

   if self._isElectromagnetic then
      self.ohmMod = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = self._basis:numBasis(),
         ghost         = {1, 1}
      }
      self.ohmMod_ptr = self.ohmMod:get(1)
      self.ohmMod_L_ptr = self.ohmMod:get(1)
      self.ohmMod_R_ptr = self.ohmMod:get(1)
      self.ohmModIdxr = self.ohmMod:genIndexer()
      self.ohmMod:clear(0.0)
   end

   -- For gyroaveraging.
   self.gyavgSlvr = tbl.gyavgSlvr
   if self.gyavgSlvr then
      self._gyavg = true
      self.phiGy  = {}
      local nmu   = self._grid:numCells(self._ndim)
      for i = 1, nmu do
         self.phiGy[i] = DataStruct.Field {
            onGrid        = self._confGrid,
            numComponents = self._confBasis:numBasis(),
            ghost         = {1, 1}
         }
      end
   end

   self._isFirst = true

   -- Timers.
   self.totalVolTime  = 0.0
   self.totalSurfTime = 0.0

   if self._positivity then
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self._grid,
         basis  = self._basis,
      }

      if self._isElectromagnetic then 
         self.fRhsVolX = DataStruct.Field {
            onGrid        = self._grid,
            numComponents = self._basis:numBasis(),
            ghost         = {1, 1}
         }
         self.fRhsVolX_ptr = self.fRhsVolX:get(1)

         self.fRhsVolX = DataStruct.Field {
            onGrid        = self._grid,
            numComponents = self._basis:numBasis(),
            ghost         = {1, 1}
         }
         self.fRhsVolX_ptr = self.fRhsVolX:get(1)

         self.fRhsSurfX = DataStruct.Field {
            onGrid        = self._grid,
            numComponents = self._basis:numBasis(),
            ghost         = {1, 1}
         }
         self.fRhsSurfX_L_ptr = self.fRhsSurfX:get(1)
         self.fRhsSurfX_R_ptr = self.fRhsSurfX:get(1)

         self.fRhsSurfX = DataStruct.Field {
            onGrid        = self._grid,
            numComponents = self._basis:numBasis(),
            ghost         = {1, 1}
         }
         self.fRhsSurfX_ptr   = self.fRhsSurfX:get(1)
         self.fRhsSurfX_L_ptr = self.fRhsSurfX:get(1)
         self.fRhsSurfX_R_ptr = self.fRhsSurfX:get(1)

         self.fRhsIdxr = self.fRhsVolX:genIndexer()
      else
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
         numComponents = self._cdim+2,
         ghost         = {1, 1},
      }
      self.positivityWeightByDirIdxr = self.positivityWeightByDir:genIndexer()
      self.positivityWeightByDir_ptr   = self.positivityWeightByDir:get(1)
      self.positivityWeightByDir_L_ptr = self.positivityWeightByDir:get(1)
      self.positivityWeightByDir_R_ptr = self.positivityWeightByDir:get(1)

      self.volTermScaleFac = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = 1,   -- Only need one number per cell.
         ghost         = {1, 1}
      }
      self.volTermScaleFacPrev = DataStruct.Field {
         onGrid        = self._grid,
         numComponents = 1,   -- Only need one number per cell.
         ghost         = {1, 1}
      }
   end
end

function Gyrokinetic:setAuxFields(auxFields)
   local potentials = auxFields[1]   -- First auxField is Field object.
   local geo        = auxFields[2]   -- Second auxField is FuncField object.

   -- Get phi.
   self.phi = potentials.phi
   if self._gyavg then 
      self.gyavgSlvr:advance(0, {self.phi}, {self.phiGy}) 
      for i=1,self._grid:numCells(self._ndim) do
         self.phiGy[i]:sync()
      end
   end

   if self._isElectromagnetic then
      -- Get electromagnetic terms.
      self.apar        = potentials.apar
      self.dApardt     = potentials.dApardt
      self.dApardtProv = auxFields[3]
   end

   -- Get magnetic geometry fields.
   self.bmag    = geo.bmag
   self.bmagInv = geo.bmagInv
   self.gradpar = geo.gradpar
   self.bdriftX = geo.bdriftX
   self.bdriftY = geo.bdriftY
   self.phiWall = geo.phiWall  -- For sheath BCs.

   if self._isFirst then
      -- Allocate pointers and indexers to field objects.

      -- Potentials.
      self.phi_ptr = self.phi:get(1)
      self.phiIdxr = self.phi:genIndexer()
      if self._isElectromagnetic then
         self.apar_ptr        = self.apar:get(1)
         self.dApardt_ptr     = self.dApardt:get(1)
         self.dApardtProv_ptr = self.dApardtProv:get(1)
         self.aparIdxr        = self.apar:genIndexer()
         self.dApardtIdxr     = self.dApardt:genIndexer()
      end

      -- For gyroaveraging.
      if self._gyavg then
         self.phiGy_ptr = {}
         self.phiGyIdxr = {}
         for i=1,self._grid:numCells(self._ndim) do
            self.phiGy_ptr[i] = self.phiGy[i]:get(1)
            self.phiGyIdxr[i] = self.phiGy[i]:genIndexer()
         end
      end

      -- Geometry.
      self.bmag_ptr    = self.bmag:get(1)
      self.bmagInv_ptr = self.bmagInv:get(1)
      self.gradpar_ptr = self.gradpar:get(1)
      self.bdriftX_ptr = self.bdriftX:get(1)
      self.bdriftY_ptr = self.bdriftY:get(1)
      self.phiWall_ptr = self.phiWall:get(1)
      self.bmagIdxr    = self.bmag:genIndexer()
      self.bmagInvIdxr = self.bmagInv:genIndexer()
      self.gradparIdxr = self.gradpar:genIndexer()
      self.bdriftXIdxr = self.bdriftX:genIndexer()
      self.bdriftYIdxr = self.bdriftY:genIndexer()
      self.phiWallIdxr = self.phiWall:genIndexer()

      self._isFirst = false   -- No longer first time.
   end
end

-- Volume integral term for use in DG scheme.
function Gyrokinetic:volTerm(w, dx, idx, f_ptr, fRhs_ptr)
   local tmStart = Time.clock()
   if self._gyavg then 
      local idmu = idx[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idx), self.phiGy_ptr[idmu])
      self.phi_ptr = self.phiGy_ptr[idmu]
   else
      self.phi:fill(self.phiIdxr(idx), self.phi_ptr)
   end
   self.bmag:fill(self.bmagIdxr(idx), self.bmag_ptr)
   self.bmagInv:fill(self.bmagInvIdxr(idx), self.bmagInv_ptr)
   self.gradpar:fill(self.gradparIdxr(idx), self.gradpar_ptr)
   self.bdriftX:fill(self.bdriftXIdxr(idx), self.bdriftX_ptr)
   self.bdriftY:fill(self.bdriftYIdxr(idx), self.bdriftY_ptr)

   local cflRate
   if self._isElectromagnetic then
      self.apar:fill(self.aparIdxr(idx), self.apar_ptr)
      if self._positivity then
         self.fRhsVolX:fill(self.fRhsIdxr(idx), self.fRhsVolX_ptr)
         self.fRhsVolX:fill(self.fRhsIdxr(idx), self.fRhsVolX_ptr)
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)

         cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), 
                             self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                             self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), self.apar_ptr:data(), 
                             f_ptr:data(), self.fRhsVolX_ptr:data(), self.fRhsVolX_ptr:data(), 
                             self.positivityWeightByDir_ptr:data())
      else
         cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), 
                             self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                             self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), self.apar_ptr:data(),
                             f_ptr:data(), fRhs_ptr:data())
      end
   else  -- Electrostatic.
      if self._positivity then 
         self.fRhsVol:fill(self.fRhsIdxr(idx), self.fRhsVol_ptr)
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)

         cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), 
                             self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                             self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                             f_ptr:data(), self.fRhsVol_ptr:data(), 
                             self.positivityWeightByDir_ptr:data())
      else
         cflRate = self._volTerm(self.charge, self.mass, w:data(), dx:data(), 
                             self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                             self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                             f_ptr:data(), fRhs_ptr:data())
      end
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return cflRate
end

-- Surface integral term for use in DG scheme
function Gyrokinetic:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, 
                              f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   local tmStart = Time.clock()
   if self._gyavg then 
      local idmu = idxr[self._ndim]
      self.phiGy[idmu]:fill(self.phiGyIdxr[idmu](idxr), self.phiGy_ptr[idmu])
      self.phi_ptr = self.phiGy_ptr[idmu]
   else
      self.phi:fill(self.phiIdxr(idxr), self.phi_ptr)
   end
   self.bmag:fill(self.bmagIdxr(idxr), self.bmag_ptr)
   self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInv_ptr)
   self.gradpar:fill(self.gradparIdxr(idxr), self.gradpar_ptr)
   self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftX_ptr)
   self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftY_ptr)

   local res
   if self._isElectromagnetic then
      self.apar:fill(self.aparIdxr(idxr), self.apar_ptr)

      if self._positivity then
         local fRhsSurf_L_ptr, fRhsSurf_R_ptr
         if dir == self._cdim + 1 then 
             self.fRhsSurfX:fill(self.fRhsIdxr(idxl), self.fRhsSurfX_L_ptr)
             self.fRhsSurfX:fill(self.fRhsIdxr(idxr), self.fRhsSurfX_R_ptr)
             fRhsSurf_L_ptr = self.fRhsSurfX_L_ptr
             fRhsSurf_R_ptr = self.fRhsSurfX_R_ptr
         else
             self.fRhsSurfX:fill(self.fRhsIdxr(idxl), self.fRhsSurfX_L_ptr)
             self.fRhsSurfX:fill(self.fRhsIdxr(idxr), self.fRhsSurfX_R_ptr)
             fRhsSurf_L_ptr = self.fRhsSurfX_L_ptr
             fRhsSurf_R_ptr = self.fRhsSurfX_R_ptr
         end
      
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxl), self.positivityWeightByDir_L_ptr)
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxr), self.positivityWeightByDir_R_ptr)

         res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), maxs,
                     self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                     self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                     self.apar_ptr:data(), self.dApardt_ptr:data(), 
                     dtApprox, self.positivityWeightByDir_L_ptr:data(), self.positivityWeightByDir_R_ptr:data(),
                     f_L_ptr:data(), f_R_ptr:data(), fRhsSurf_L_ptr:data(), fRhsSurf_R_ptr:data())
      else
         self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProv_ptr)
         self.ohmMod:fill(self.ohmModIdxr(idxl), self.ohmMod_L_ptr)
         self.ohmMod:fill(self.ohmModIdxr(idxr), self.ohmMod_R_ptr)
         res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), maxs,
                     self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                     self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                     self.apar_ptr:data(), self.dApardt_ptr:data(), self.dApardtProv_ptr:data(), 
                     f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data(), 
                     self.ohmMod_L_ptr:data(), self.ohmMod_R_ptr:data())
      end
   else   -- Electrostatic.
      if self._positivity then
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxl), self.positivityWeightByDir_L_ptr)
         self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idxr), self.positivityWeightByDir_R_ptr)

         self.fRhsSurf:fill(self.fRhsIdxr(idxl), self.fRhsSurf_L_ptr)
         self.fRhsSurf:fill(self.fRhsIdxr(idxr), self.fRhsSurf_R_ptr)

         res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), maxs,
                     self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                     self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                     dtApprox, self.positivityWeightByDir_L_ptr:data(), self.positivityWeightByDir_R_ptr:data(),
                     f_L_ptr:data(), f_R_ptr:data(), self.fRhsSurf_L_ptr:data(), self.fRhsSurf_R_ptr:data())
      else
         res = self._surfTerms[dir](self.charge, self.mass, wr:data(), dxr:data(), maxs, 
                     self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                     self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                     f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data())
      end
   end
   self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   return res
end

function Gyrokinetic:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f_ptr, fRefl_ptr)
   self.phi:fill(self.phiIdxr(idx), self.phi_ptr)
   self.phiWall:fill(self.phiWallIdxr(idx), self.phiWall_ptr)
   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, 
                                            self.phi_ptr:data(), self.phiWall_ptr:data(), f_ptr:data(), fRefl_ptr:data())
end

-- Step2 volume integral term for use in DG scheme (EM only).
function Gyrokinetic:volTermStep2(w, dx, idx, f_ptr, fRhs_ptr)
   local tmStart = Time.clock()
   self.dApardt:fill(self.dApardtIdxr(idx), self.dApardt_ptr)
   local cflRate
   if self._positivity then
      self.fRhsVolX:fill(self.fRhsIdxr(idx), self.fRhsVolX_ptr)
      self.fRhsSurfX:fill(self.fRhsIdxr(idx), self.fRhsSurfX_ptr)
      self.positivityWeightByDir:fill(self.positivityWeightByDirIdxr(idx), self.positivityWeightByDir_ptr)

      cflRate = self._volTermStep2(self.charge, self.mass, w:data(), dx:data(), 
                                   self.dApardt_ptr:data(), 
                                   f_ptr:data(), self.fRhsVolX_ptr:data(),
                                   self.positivityWeightByDir_ptr:data())
   else
      cflRate = self._volTermStep2(self.charge, self.mass, w:data(), dx:data(), 
                                   self.dApardt_ptr:data(), 
                                   f_ptr:data(), fRhs_ptr:data())
   end
   self.totalVolTime = self.totalVolTime + (Time.clock()-tmStart)
   return cflRate
end

-- Step2 surface integral term for use in DG scheme (EM only, vpar dir only).
function Gyrokinetic:surfTermStep2(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, 
                                   f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   local res
   if self._positivity or self._polyOrder > 1 then
      -- no separate step 2 surfTerm kernel, just do entire vpar surf term in step 2
      return self:surfTerm(dir, dtApprox, wl, wr, dxl, dxr, maxs, idxl, idxr, 
                                      f_L_ptr, f_R_ptr, fRhs_L_ptr, fRhs_R_ptr)
   else
      local tmStart = Time.clock()
      self.phi:fill(self.phiIdxr(idxr), self.phi_ptr)
      self.bmag:fill(self.bmagIdxr(idxr), self.bmag_ptr)
      self.bmagInv:fill(self.bmagInvIdxr(idxr), self.bmagInv_ptr)
      self.gradpar:fill(self.gradparIdxr(idxr), self.gradpar_ptr)
      self.bdriftX:fill(self.bdriftXIdxr(idxr), self.bdriftX_ptr)
      self.bdriftY:fill(self.bdriftYIdxr(idxr), self.bdriftY_ptr)
      self.apar:fill(self.aparIdxr(idxr), self.apar_ptr)
      self.dApardt:fill(self.dApardtIdxr(idxr), self.dApardt_ptr)
      self.dApardtProv:fill(self.dApardtIdxr(idxr), self.dApardtProv_ptr)

      res = self._surfTermsStep2(self.charge, self.mass, wr:data(), dxr:data(), maxs,
                                      self.bmag_ptr:data(), self.bmagInv_ptr:data(), self.gradpar_ptr:data(), 
                                      self.bdriftX_ptr:data(), self.bdriftY_ptr:data(), self.phi_ptr:data(), 
                                      self.apar_ptr:data(), self.dApardt_ptr:data(), self.dApardtProv_ptr:data(), 
                                      f_L_ptr:data(), f_R_ptr:data(), fRhs_L_ptr:data(), fRhs_R_ptr:data())
      self.totalSurfTime = self.totalSurfTime + (Time.clock()-tmStart)
   end

   return res
end

function Gyrokinetic:clearRhsTerms()
   if self._positivity then
      if self._isElectromagnetic then
         self.fRhsVolX:clear(0.0)
         self.fRhsVolX:clear(0.0)
         self.fRhsSurfX:clear(0.0)
         self.fRhsSurfX:clear(0.0)
      else    -- Electrostatic.
         self.fRhsVol:clear(0.0)
         self.fRhsSurf:clear(0.0)
      end
   end
end

-- When using positivity algorithm, different parts of RHS are stored separately.
-- here we combine the parts, with some rescaling of the volume term.
function Gyrokinetic:getPositivityRhs(tCurr, dtApprox, fIn, fRhs)
   if self._isElectromagnetic then
      weightDirs = {}
      for d = 1, self._cdim do -- configuration space dirs
         weightDirs[d] = d
      end
      --if dtApprox > 0 then self.posRescaler:rescaleVolTerm(tCurr, dtApprox*1.05, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurfX, self.fRhsVolX) end

      -- accumulate Rhs terms into total fRhs. 
      -- note that fRhsVolX is only the partial vpar volume term, the part that doesn't depend on dA/dt. this is because this term is needed for Ohm's law.
      -- before the step 2 update we will remove this partial fRhsVolX from fRhs. 
      -- then the step 2 update will accumulate the dA/dt vol term into fRhsVolX, and we will add the entire fRhsVolX into fRhs (in getPositivityRhsStep2).
      fRhs:accumulate(1.0, self.fRhsVolX)
      fRhs:scaleByCell(self.volTermScaleFac)
      fRhs:accumulate(1.0, self.fRhsSurfX)
   else
      weightDirs = {}
      for d = 1, self._cdim+1 do -- configuration space dirs + vpar dir
         weightDirs[d] = d
      end
      -- fIn + fac*dt*fVol + dt*fSurf > 0.
      -- Rescale volume term by fac, and add to surface term fRhs = fRhsSurf.
      if dtApprox > 0 then self.posRescaler:rescaleVolTerm(tCurr, dtApprox*1.05, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurf, self.fRhsVol) end

      fRhs:accumulate(1.0, self.fRhsSurf, 1.0, self.fRhsVol)
   end 
end

function Gyrokinetic:clearRhsTermsStep2(fRhs)
   if self._positivity and self._isElectromagnetic then
      -- fRhsVolX that was accumulated into fRhs was only partial, because it was needed for Ohm's law. 
      -- now subtract it, so that later we can accumulate the entire fRhsVolX into fRhs.
      --fRhs:accumulate(-1.0, self.fRhsVolX)
      self.fRhsVolX:clear(0.0)
   end
end

function Gyrokinetic:getPositivityRhsStep2(tCurr, dtApprox, fIn, fRhs)
   if self._isElectromagnetic then
      weightDirs = {self._cdim + 1} -- only vpar direction
      --if dtApprox > 0 then self.posRescaler:rescaleVolTerm(tCurr, dtApprox*1.05, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurfX, self.fRhsVolX) end
      
      self.volTermScaleFacPrev:copy(self.volTermScaleFac)
      if dtApprox > 0 then 
         self.posRescaler:calcVolTermRescale(tCurr, dtApprox*1.05, fIn, self.positivityWeightByDir, weightDirs, self.fRhsSurfX, self.fRhsVolX, self.volTermScaleFac) 
      end

      self.fRhsVolX:scaleByCell(self.volTermScaleFacPrev)
      fRhs:combine(1.0, self.fRhsSurfX, 1.0, self.fRhsVolX)
   end 
end

function Gyrokinetic:setPositivityWeights(cflRateByCell)
   -- set total weight = positivityWeightByDir[0] = cflRateByCell in each cell
   -- other elements will be overwritten in kernels
   self.positivityWeightByDir:clear(1.0)
   self.positivityWeightByDir:scaleByCell(cflRateByCell)
end

function Gyrokinetic:sync()
   if self._positivity then 
      self.positivityWeightByDir:sync()
   end
end

function Gyrokinetic:read(name)
   self.positivityWeightByDir:read(
      string.format("%s_cflRateByCell_restart.bp", name), false)
   self.positivityWeightByDir:sync()
end

function Gyrokinetic:write(name, tm, frame)
   self.positivityWeightByDir:write(
      string.format("%s_cflRateByCell_restart.bp", name), tm, frame, false)
end

return {GkEq = Gyrokinetic}
