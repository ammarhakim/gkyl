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
   -- setup base object
   BraginskiiViscosityDiffusion.super.init(self, tbl)

   local pfx = "Updater.BraginskiiViscosityDiffusion: "

   self._onGrid = assert(tbl.onGrid, pfx.." Must provide 'onGrid'.")

   self._gasGamma = assert(tbl.gasGamma, pfx .. "Must provide 'gasGamma'.")

   self._nFluids = assert(tbl.numFluids, pfx .. "Must provide 'numFluids'.")

   self._eta = assert(tbl.eta, pfx.."Must provide 'eta'")

   self._hasHeating = tbl.hasHeating ~= nil and tbl.hasHeating or false

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.",
                        pfx, tbl._coordinate))

   assert(self._gasGamma==5./3., pfx .. " gasGamma must be 5/3.")
end

function BraginskiiViscosityDiffusion:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local nFluids = self._nFluids

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   -- Indices.
   local ndim = grid:ndim()
   local idxp = Lin.IntVec(grid:ndim())
   local idxm = Lin.IntVec(grid:ndim())

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

               local eta = self._eta
               local etaP = self._eta
               local etaM = self._eta

               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               -- Compute momentum diffusion as grad(eta * grad(V)).
               local dx = grid:dx(d)
               for c=2,4 do
                  bufPtr[c] = (etaP*(fldPtrP[c]-fldPtr[c]) -
                               etaM*(fldPtr[c]-fldPtrM[c])) / (dx*dx)
               end

               -- Compute viscous heating eta*grad(V):grad(V) as
               -- eta*|grad(V)|^2.
               if self._hasHeating then
                  local eta = self._eta
                  for c=2,4 do
                     qVis = qVis + eta * ( (fldPtrP[c]-fldPtrM[c])/(2*dx) )^2
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

            -- Radial diffusion.
            if true then
               local d = 1

               local dr = grid:dx(1)
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
               fld:fill(fldIdxr(idx), fldPtr)
               fld:fill(fldIdxr(idxp), fldPtrP)
               fld:fill(fldIdxr(idxm), fldPtrM)
               local eta = self._eta
               local etaP = self._eta
               local etaM = self._eta

               grid:setIndex(idx)
               grid:cellCenter(xc)
               grid:setIndex(idxp)
               grid:cellCenter(xp)
               grid:setIndex(idxm)
               grid:cellCenter(xm)
               local r = xc[1]
               local rp = xp[1]
               local rm = xm[1]

               local rpH = 0.5 * (r+rp)
               local rmH = 0.5 * (r+rm)
               -- FIXME Is it better to compute eta with average B, etc.?
               local etaMH = 0.5 * (eta+etaM)
               local etaPH = 0.5 * (eta+etaP)

               bufPtr[2] =
                  (etaPH*rpH*(fldPtrP[2]-fldPtr [2]) -
                   etaMH*rmH*(fldPtr [2]-fldPtrM[2])) / (dr*dr*r) -
                  eta*fldPtr[2]/(r*r)
               bufPtr[3] =
                  (etaPH*(rpH^3)*(fldPtrP[3]/rp-fldPtr [3]/r) -
                   etaMH*(rmH^3)*(fldPtr [3]/r -fldPtrM[3]/rm)) / (dr*dr*r)
               bufPtr[4] =
                  (etaPH*rpH*(fldPtrP[4]-fldPtr [4]) -
                   etaMH*rmH*(fldPtr [4]-fldPtrM[4])) / (dr*dr*r)

               if self._hasHeating then
                  local eta = self._eta
                  qVis = qVis + eta * (fldPtrP[3]/rp-fldPtrM[3]/rm) * 0.5 / dr
               end
            end

            -- Axial diffusion.
            if true then
               d=3

               local dz = grid:dx(3)
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
               fld:fill(fldIdxr(idx), fldPtr)
               fld:fill(fldIdxr(idxp), fldPtrP)
               fld:fill(fldIdxr(idxm), fldPtrM)
               local eta = self._eta
               local etaP = self._eta
               local etaM = self._eta

               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               for c=2,4 do
                  bufPtr[c] = bufPtr[c] +
                              (etaPH*(fldPtrP[c]-fldPtr [c]) -
                               etaMH*(fldPtr[c] -fldPtrM[c])) / (dz*dz)
               end

               if self._hasHeating then
                  grid:setIndex(idx)
                  grid:cellCenter(xc)
                  qVis = qVis + eta * (fldPtrP[3]-fldPtrM[3]) / r * 0.5 / dz
               end
            end

            bufPtr[5] = qVis
         end
      else
         assert(false)
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

         if self._hasHeating then
            local keOld = 0.5*(fldPtr[2]^2+fldPtr[3]^2+fldPtr[4]^2) / fldPtr[1]
            for c=2,4 do
               fldPtr[c] = fldPtr[c] + bufPtr[c]
            end
            local keNew = 0.5*(fldPtr[2]^2+fldPtr[3]^2+fldPtr[4]^2) / fldPtr[1]
            fldPtr[5] = fldPtr[5]+keNew-keOld+bufPtr[5]
         else
            for c=2,4 do fldPtr[c] = fldPtr[c] + bufPtr[c] end
         end
      end
   end

   return status, dtSuggested
end

function BraginskiiViscosityDiffusion:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiViscosityDiffusion
