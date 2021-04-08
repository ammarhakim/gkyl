-- Gkyl ------------------------------------------------------------------------
--
-- Updater to apply isotropic, diffusion-like viscosity instead of full
-- Braginskii viscosity tensor.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local BraginskiiViscosityDiffusion = Proto(UpdaterBase)

function BraginskiiViscosityDiffusion:init(tbl)
   BraginskiiViscosityDiffusion.super.init(self, tbl)

   local pfx = "Updater.BraginskiiViscosityDiffusion: "

   self._onGrid = assert(tbl.onGrid, pfx.." Must provide 'onGrid'.")
   self._cfl = tbl.cfl and tbl.cfl or 1

   self._nFluids = assert(tbl.numFluids, pfx .. "Must provide 'numFluids'.")

   self._eta = assert(tbl.eta, pfx.."Must provide 'eta'")
   assert(#self._eta==self._nFluids,
          pfx.."Number of entries in 'eta' must match numFluids.")

   self._hasHeating = tbl.hasHeating ~= nil and tbl.hasHeating or false

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.",
                        pfx, tbl._coordinate))
end

-- TODO Nicer, more accurate time-step size calculation.
local suggestDt = function(eta, ndim, dx, cfl)
      local eta__dx2_sum = 0
      for d = 1, ndim do
         eta__dx2_sum = eta__dx2_sum + eta / ( (dx[d])^2 )
      end
      return cfl * 0.5 / eta__dx2_sum
end

function BraginskiiViscosityDiffusion:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local nFluids = self._nFluids

   -- Indices.
   local ndim = grid:ndim()
   local idxp = Lin.IntVec(grid:ndim())
   local idxm = Lin.IntVec(grid:ndim())

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   local dx = {grid:dx(1), grid:dx(2), grid:dx(3)}
   for s = 1, nFluids do
      dtSuggested = suggestDt(self._eta[s], ndim, dx, self._cfl)
      status = dt <= dtSuggested
      if not status then return false, dtSuggested end
   end

   -- Comptue grad_para(T) ain internal cells.
   for s = 1, nFluids do
      local fld = outFld[s]
      local fldIdxr = fld:genIndexer()
      local fldPtr = fld:get(1)
      local fldPtrP = fld:get(1)
      local fldPtrM = fld:get(1)

      local buf = inFld[s]
      local bufIdxr = buf:genIndexer()
      local bufPtr = buf:get(1)
      local bufPtrP = buf:get(1)
      local bufPtrM = buf:get(1)

      local localRange = fld:localRange()

      -- Compute rhs.
      if self._coordinate=="cartesian" then
         for idx in localRange:rowMajorIter() do
            local qVis = 0
            for d = 1, ndim do
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
              
               fld:fill(fldIdxr(idx), fldPtr)
               fld:fill(fldIdxr(idxp), fldPtrP)
               fld:fill(fldIdxr(idxm), fldPtrM)
               buf:fill(bufIdxr(idx), bufPtr)

               local eta = self._eta[s]
               local etaP = self._eta[s]
               local etaM = self._eta[s]
               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               local rho = fldPtrP[1]
               local rhoP = fldPtrP[1]
               local rhoM = fldPtrM[1]
               local rhoPH = 0.5 * (rho+rhoP)
               local rhoMH = 0.5 * (rho+rhoM)

               -- Compute momentum diffusion as div(eta * rho * grad(V)).
               local dx = grid:dx(d)
               for c=2,4 do
                  bufPtr[c] =
                     ( rhoPH*etaPH*(fldPtrP[c]/rhoP - fldPtr [c]/rho ) -
                       rhoMH*etaMH*(fldPtr [c]/rho  - fldPtrM[c]/rhoM) ) /
                     (dx*dx)
               end

               -- Compute viscous heating eta * rho * grad(V) : grad(V) as
               -- eta * rho * |grad(V)|^2.
               if self._hasHeating then
                  for c=2,4 do
                     qVis = qVis + eta * rho * 
                            ( (fldPtrP[c]/rhoP-fldPtrM[c]/rhoM)/(2*dx) )^2
                  end
               end
            end
            bufPtr[5] = qVis
         end
      elseif self._coordinate=="axisymmetric" then
         local xc = Lin.Vec(ndim)
         local xp = Lin.Vec(ndim)
         local xm = Lin.Vec(ndim)
         for idx in localRange:rowMajorIter() do
            local qVis = 0

            fld:fill(fldIdxr(idx), fldPtr)
            local eta = self._eta[s]
            local rho = fldPtrP[1]
            grid:setIndex(idx)
            grid:cellCenter(xc)
            local r = xc[1]

            -- Radial diffusion.
            if true then
               local d = 1

               local dr = grid:dx(1)
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
               fld:fill(fldIdxr(idxp), fldPtrP)
               fld:fill(fldIdxr(idxm), fldPtrM)

               local etaP = self._eta[s]
               local etaM = self._eta[s]
               local etaMH = 0.5 * (eta+etaM)
               local etaPH = 0.5 * (eta+etaP)

               grid:setIndex(idxp)
               grid:cellCenter(xp)
               local rp = xp[1]
               grid:setIndex(idxm)
               grid:cellCenter(xm)
               local rm = xm[1]
               local rpH = 0.5 * (r+rp)
               local rmH = 0.5 * (r+rm)

               local rhoP = fldPtrP[1]
               local rhoM = fldPtrM[1]
               local rhoPH = 0.5 * (rho+rhoP)
               local rhoMH = 0.5 * (rho+rhoM)

               bufPtr[2] =
                  ( rhoPH*etaPH*rpH * (fldPtrP[2]/rhoP - fldPtr [2]/rho ) -
                    rhoMH*etaMH*rmH * (fldPtr [2]/rho  - fldPtrM[2]/rhoM) ) /
                  (dr*dr*r) -
                  eta*fldPtr[2]/(r*r)
               bufPtr[3] =
                  (rhoPH*etaPH*(rpH^2)*(fldPtrP[3]/rhoP/rp-fldPtr [3]/rho /r ) -
                   rhoMH*etaMH*(rmH^2)*(fldPtr [3]/rho /r -fldPtrM[3]/rhoM/rm))/
                  (dr*dr*r)
               bufPtr[4] =
                  ( rhoPH*etaPH*rpH * (fldPtrP[4]/rhoP - fldPtr [4]/rho ) -
                    rhoMH*etaMH*rmH * (fldPtr [4]/rho  - fldPtrM[4]/rhoM) ) /
                  (dr*dr*r)

               if self._hasHeating then
                  local eta = self._eta[s]
                  qVis = qVis +
                         rho * eta *
                         ( fldPtrP[3]/rhoP/rp-fldPtrM[3]/rhoM/rm ) / (2*dr)
               end
            end

            -- Axial diffusion.
            if ndim>=3 then
               d=3

               local dz = grid:dx(3)
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
               fld:fill(fldIdxr(idxp), fldPtrP)
               fld:fill(fldIdxr(idxm), fldPtrM)

               local etaP = self._eta[s]
               local etaM = self._eta[s]
               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               local rhoP = fldPtrP[1]
               local rhoM = fldPtrM[1]
               local rhoPH = 0.5 * (rho+rhoP)
               local rhoMH = 0.5 * (rho+rhoM)

               for c=2,4 do
                  bufPtr[c] = bufPtr[c] +
                              (rhoPH*etaPH*(fldPtrP[c]/rhoP-fldPtr [c]/rho ) -
                               rhoMH*etaMH*(fldPtr[c] /rho -fldPtrM[c]/rhoM)) /
                              (dz*dz)
               end

               if self._hasHeating then
                  grid:setIndex(idx)
                  grid:cellCenter(xc)
                  qVis = qVis +
                         rho*eta* (fldPtrP[3]/rhoP-fldPtrM[3]/rhoM) / r / (2*dz)
               end
            end

            bufPtr[5] = qVis
         end
      else
         assert(false, "Coordinate type "..self._coordinate.." not supported.")
      end
   end

   -- Apply rhs.
   for s = 1, nFluids do
      local fld = outFld[s]
      local fldIdxr = fld:genIndexer()
      local fldPtr = fld:get(1)

      local buf = inFld[s]
      local bufIdxr = buf:genIndexer()
      local bufPtr = buf:get(1)

      local localRange = fld:localRange()

      for idx in localRange:rowMajorIter() do
         fld:fill(fldIdxr(idx), fldPtr)
         buf:fill(bufIdxr(idx), bufPtr)

         local keOld = 0.5*(fldPtr[2]^2+fldPtr[3]^2+fldPtr[4]^2) / fldPtr[1]

         for c=2,4 do fldPtr[c] = fldPtr[c] + bufPtr[c] * dt end

         local keNew = 0.5*(fldPtr[2]^2+fldPtr[3]^2+fldPtr[4]^2) / fldPtr[1]

         local dEnergy = keNew - keOld
         if self._hasHeating then
            dEnergy = dEnergy + bufPtr[5] * dt
         end

         fldPtr[5] = fldPtr[5] + dEnergy
      end
   end

   return status, dtSuggested
end

function BraginskiiViscosityDiffusion:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiViscosityDiffusion
