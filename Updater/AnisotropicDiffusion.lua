-- Gkyl: Anisotropic diffusion of a scalar parallel and perpendicular to the
-- background magnetic field:
--
-- dT/dt = -div(q), where q = -kappaPara*gradPara(T)-kappaPerp*gradPerp(T).
--
-- Mainly for testing purpose. The verified schemes are used to construct more
-- complicated heat conduction, viscosity updaters, etc.
--
-- An explicit two-step algorithm is used:
--   First, calculate grad(T) and then q;
--   Second, calculate div(q), and add it onto T.
--
-- Major TODOs:
-- - More accurate time-step size estimation.
-- - Limiters following Sharma and Hammett 2007 JCP.
--
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local AnisotropicDiffusion = Proto(UpdaterBase)

function AnisotropicDiffusion:init(tbl)
   AnisotropicDiffusion.super.init(self, tbl)
   local pfx = "Updater.AnisotropicDiffusion: "

   self._onGrid = assert(tbl.onGrid, pfx.."Must provide 'onGrid'.")
   self._cfl = tbl.cfl and tbl.cfl or 1

   self._kappaPara = assert(tbl.kappaPara, pfx.."Must provide 'kappaPara'")
   self._kappaPerp = assert(tbl.kappaPerp, pfx.."Must provide 'kappaPerp'")

   self._component = tbl.component and tbl.component or 1

   self._timeStepper = tbl.timeStepper~=nil and tbl.timeStepper or
                       "symmetric-node-center"
   assert(self._timeStepper=="symmetric-cell-center" or
          self._timeStepper=="symmetric-node-center",
          pfx.."timeStepper '"..self._timeStepper.."' is not supported.")
end

local function isNaN( v ) return type( v ) == "number" and v ~= v end

function AnisotropicDiffusion:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt

   local ndim = grid:ndim()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())

   local temp = outFld[1]
   local tempIdxr = temp:genIndexer()
   local tempPtr = temp:get(1)
   local tempPtrP = temp:get(1)
   local tempPtrM = temp:get(1)

   local emf = outFld[2]
   local emfIdxr = emf:genIndexer()
   local emfPtr = emf:get(1)

   local buf = inFld[1]
   local bufIdxr = buf:genIndexer()
   local bufPtr = buf:get(1)
   local bufPtrP = buf:get(1)
   local bufPtrM = buf:get(1)

   local localRange = temp:localRange()

   -- Check time-step size.
   local kappaPara = self._kappaPara
   local kappaPerp = self._kappaPerp

   -- TODO Nicer, more accurate time-step size calculation.
   local kappaMax = math.max(kappaPara, kappaPerp)
   local kappa__dx2_sum = 0
   for d = 1, ndim do
      kappa__dx2_sum = kappa__dx2_sum + kappaMax / ((grid:dx(d))^2)
   end

   local dtSuggested = self._cfl * 0.5 / kappa__dx2_sum
   local status = dt <= dtSuggested
   if not status then
      return false, dtSuggested
   end

   local c = self._component

   if self._timeStepper=="symmetric-cell-center" then

      -- Comptue grad(T) in internal + one ghost cell centers.
      local localExt1Range = localRange:extend(1, 1)
      for idx in localExt1Range:rowMajorIter() do
         for d = 1, ndim do
            idx:copyInto(idxp)
            idx:copyInto(idxm)
            idxp[d] = idx[d]+1
            idxm[d] = idx[d]-1

            buf:fill(bufIdxr(idx), bufPtr)
            temp:fill(tempIdxr(idxp), tempPtrP)
            temp:fill(tempIdxr(idxm), tempPtrM)

            bufPtr[d] = (tempPtrP[c] - tempPtrM[c]) * 0.5 / grid:dx(d)
         end
      end

      -- Compute q = q_para + q_perp at cell centers.
      -- q_para = -kappa_para*gradPara(T), q_perp = -kappa_perp*gradPerp(T).
      for idx in localExt1Range:rowMajorIter() do
         emf:fill(emfIdxr(idx), emfPtr)
         buf:fill(bufIdxr(idx), bufPtr)

         local bx = emfPtr[4]
         local by = emfPtr[5]
         local bz = emfPtr[6]
         local bmag = math.sqrt(bx*bx + by*by + bz*bz)
         bx = bx / bmag
         by = by / bmag
         bz = bz / bmag
         assert(bmag>0, "Zero B field detected!")

         local bDotGradT = bx*bufPtr[1]
         if ndim>1 then bDotGradT = bDotGradT + by*bufPtr[2] end
         if ndim>2 then bDotGradT = bDotGradT + bz*bufPtr[3] end

         local gradParaTx = bx * bDotGradT
         gradPerpTx = bufPtr[1] - gradParaTx
         bufPtr[1] = - kappaPara*gradParaTx - kappaPerp*gradPerpTx

         if ndim>1 then
            local gradParaTy = by * bDotGradT
            gradPerpTy = bufPtr[2] - gradParaTy
            bufPtr[2] = - kappaPara*gradParaTy - kappaPerp*gradPerpTy
         end

         if ndim>2 then
            local gradParaTz = bz * bDotGradT
            gradPerpTz = bufPtr[3] - gradParaTz
            bufPtr[3] = - kappaPara*gradParaTz - kappaPerp*gradPerpTz
         end
      end

      for idx in localRange:rowMajorIter() do
         local divq = 0
         for d = 1, ndim do
            idx:copyInto(idxp)
            idx:copyInto(idxm)
            idxp[d] = idx[d]+1
            idxm[d] = idx[d]-1

            buf:fill(bufIdxr(idxp), bufPtrP)
            buf:fill(bufIdxr(idxm), bufPtrM)

            divq = divq + (bufPtrP[d] - bufPtrM[d]) * 0.5 / grid:dx(d)
         end

         buf:fill(bufIdxr(idx), bufPtr)
         bufPtr[4] = divq
      end

   elseif self._timeStepper=="symmetric-node-center" then

      -- Compute grad(T) on nodes (cell-corners).
      -- The i-th node here is defined as the lower corner of the i-th cell,
      -- therefore the node's adjacent cells have cell indices i-1 and i.
      local localNodeRange = localRange:extend(0,1)
      for idx in localNodeRange:rowMajorIter() do
         buf:fill(bufIdxr(idx), bufPtr)

         if ndim==1 then
            local d = 1
            idx:copyInto(idxp)
            idx:copyInto(idxm)
            idxp[d] = idx[d]
            idxm[d] = idx[d]-1
            temp:fill(tempIdxr(idxp), tempPtrP)
            temp:fill(tempIdxr(idxm), tempPtrM)
            bufPtr[d] = (tempPtrP[c] - tempPtrM[c]) / grid:dx(d)
         elseif ndim==2 then
            local subIterDirs = {{2}, {1}}
            for d=1,ndim do  -- Gradient direction.
               local gradT = 0  -- Grad(T) along this direction.

               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]
               idxm[d] = idx[d]-1

               -- Add contributions from cell-center values of cells sharing
               -- this corner.
               local d1 = subIterDirs[d][1]
               for _,commonShift1 in ipairs({-1, 0}) do
                  idxp[d1] = idx[d1] + commonShift1
                  idxm[d1] = idx[d1] + commonShift1
                  temp:fill(tempIdxr(idxp), tempPtrP)
                  temp:fill(tempIdxr(idxm), tempPtrM)
                  gradT = gradT + (tempPtrP[c] - tempPtrM[c])
               end
               bufPtr[d] = gradT * 0.5 / grid:dx(d) 
            end
         elseif ndim==3 then
            local subIterDirs = {{2,3}, {1,3}, {1,2}}
            for d=1,ndim do  -- Gradient direction.
               local gradT = 0  -- Grad(T) along this direction.

               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]
               idxm[d] = idx[d]-1

               -- Add contributions from cell-center values of cells sharing
               -- this corner.
               local d1 = subIterDirs[d][1]
               local d2 = subIterDirs[d][2]
               for _,commonShift1 in ipairs({-1, 0}) do
                  for _,commonShift2 in ipairs({-1, 0}) do
                     idxp[d1] = idx[d1] + commonShift1
                     idxm[d1] = idx[d1] + commonShift1
                     idxp[d2] = idx[d2] + commonShift2
                     idxm[d2] = idx[d2] + commonShift2
                     temp:fill(tempIdxr(idxp), tempPtrP)
                     temp:fill(tempIdxr(idxm), tempPtrM)
                     gradT = gradT + (tempPtrP[c] - tempPtrM[c])
                  end
               end
               bufPtr[d] = gradT * 0.25 / grid:dx(d) 
            end -- ndim==3 ends.
         end -- Loop over ndim==1,2,3 ends.
      end -- Node-center grad(T) computation ends.

      -- Compute q on nodes (cell corners).
      for idx in localNodeRange:rowMajorIter() do
         -- Compute B field at cell corners.
         -- The value at the i-th node is defined as an average of values at
         -- centers of all cells sharing this node.
         local bx, by, bz = 0, 0, 0

         local nPts = 0
         local xshifts = {-1, 0}
         local yshifts = ndim>1 and {-1, 0} or {-1}
         local zshifts = ndim>2 and {-1, 0} or {-1}
         for _,zshift in ipairs(zshifts) do
            for _,yshift in ipairs(yshifts) do
               for _,xshift in ipairs(xshifts) do
                  idx:copyInto(idxp)
                  if ndim>2 then idxp[3] = idx[3]+zshift end
                  if ndim>1 then idxp[2] = idx[2]+yshift end
                  idxp[1] = idx[1]+xshift

                  emf:fill(emfIdxr(idxp), emfPtr)
                  bx = bx + emfPtr[4]
                  by = by + emfPtr[5]
                  bz = bz + emfPtr[6]

                  nPts = nPts + 1
               end
            end
         end
         bx = bx / nPts
         by = by / nPts
         bz = bz / nPts

         local bmag = math.sqrt(bx*bx + by*by + bz*bz)
         bx = bx / bmag
         by = by / bmag
         bz = bz / bmag
         assert(bmag>0, "Zero B field detected!")

         -- Compute gradParaT and gradPerpT at cell corners.
         buf:fill(bufIdxr(idx), bufPtr)
         local bDotGradT = bx*bufPtr[1] + by*bufPtr[2] + bz*bufPtr[3]

         local gradParaTx = bx * bDotGradT
         local gradParaTy = by * bDotGradT
         local gradParaTz = bz * bDotGradT

         gradPerpTx = bufPtr[1] - gradParaTx
         gradPerpTy = bufPtr[2] - gradParaTy
         gradPerpTz = bufPtr[3] - gradParaTz

         buf:fill(bufIdxr(idx), bufPtr)

         bufPtr[1] = - kappaPara*gradParaTx - kappaPerp*gradPerpTx
         bufPtr[2] = - kappaPara*gradParaTy - kappaPerp*gradPerpTy
         bufPtr[3] = - kappaPara*gradParaTz - kappaPerp*gradPerpTz
      end  -- Node-center q computation ends.

      -- Compute div(q) at cell centers.
      for idx in localRange:rowMajorIter() do
         buf:fill(bufIdxr(idx), bufPtr)
         local divq = 0

         if ndim==1 then
            local d = 1
            idx:copyInto(idxp)
            idx:copyInto(idxm)
            idxp[d] = idx[d]+1
            idxm[d] = idx[d]
            buf:fill(bufIdxr(idxp), bufPtrP)
            buf:fill(bufIdxr(idxm), bufPtrM)
            divq = divq + (bufPtrP[d]-bufPtrM[d]) / grid:dx(d)
         elseif ndim==2 then
            local subIterDirs = {{2}, {1}}
            for d=1,ndim do  -- Gradient direction.
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]

               -- Add contributions from node-center values.
               local d1 = subIterDirs[d][1]
               for _,commonShift1 in ipairs({0, 1}) do
                  idxp[d1] = idx[d1] + commonShift1
                  idxm[d1] = idx[d1] + commonShift1
                  buf:fill(bufIdxr(idxp), bufPtrP)
                  buf:fill(bufIdxr(idxm), bufPtrM)
                  divq = divq + (bufPtrP[d]-bufPtrM[d])*0.5/grid:dx(d)
               end
            end
         elseif ndim==3 then
            local subIterDirs = {{2,3}, {1,3}, {1,2}}
            for d=1,ndim do  -- Gradient direction.
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]

               -- Add contributions from node-center values.
               local d1 = subIterDirs[d][1]
               local d2 = subIterDirs[d][2]
               for _,commonShift1 in ipairs({0, 1}) do
                  for _,commonShift2 in ipairs({0, 1}) do
                     idxp[d1] = idx[d1] + commonShift1
                     idxm[d1] = idx[d1] + commonShift1
                     idxp[d2] = idx[d2] + commonShift2
                     idxm[d2] = idx[d2] + commonShift2
                     buf:fill(bufIdxr(idxp), bufPtrP)
                     buf:fill(bufIdxr(idxm), bufPtrM)
                     divq = divq + (bufPtrP[d]-bufPtrM[d])*0.25/grid:dx(d)
                  end
               end
            end -- Loop over gradient directions.
         end -- Loop over ndim==1,2,3 ends.
         bufPtr[4] = divq 
      end -- div(q) computation ends.

   end  -- timeStepper handling in divq calculation ends.

   -- Add div(q) to energy.
   for idx in localRange:rowMajorIter() do
      temp:fill(tempIdxr(idx), tempPtr)
      buf:fill(bufIdxr(idx), bufPtr)
      tempPtr[1] = tempPtr[1] - dt * bufPtr[4]
   end

   return status, dtSuggested
end

function AnisotropicDiffusion:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return AnisotropicDiffusion
