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
-- - Non-uniform grid support.
--
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local xsys = require "xsys"

local AnisotropicDiffusion = Proto(UpdaterBase)

function AnisotropicDiffusion:init(tbl)
   AnisotropicDiffusion.super.init(self, tbl)
   local pfx = "Updater.AnisotropicDiffusion: "
   self._pfx = pfx

   self._onGrid = assert(tbl.onGrid, pfx.."Must provide 'onGrid'.")
   self._cfl = tbl.cfl and tbl.cfl or 1

   self._kappaMode = tbl.kappaMode and tbl.kappaMode or "constant"
   if self._kappaMode=="constant" then
      local pfxx = pfx.."For kappaMode=='constant', "
      self._kappaPara = assert(tbl.kappaPara, pfxx.."Must provide 'kappaPara'.")
      self._kappaPerp = assert(tbl.kappaPerp, pfxx.."Must provide 'kappaPerp'.")
   elseif self._kappaMode=="field" then
      self._kappaField = assert(self.kappaField,
                                pfx.."For kappaMode=='field', "..
                                "Must provide 'kappaField'.")
   elseif self._kappaMode=="function" then
      self._kappaFunction = assert(tbl.kappaFunction,
                                   pfx.."For kappaMode=='function', "..
                                   "Must provide 'kappaFunction'.")
   end

   -- Index in the input that contains the variable whose diffusion is to be
   -- computed.
   self._components = tbl.components and tbl.components or {1}
   -- Indices in the buffer to store temporary grad(T) vector which will then be
   -- overwritten by q vector.
   self._componentsBufQ = tbl.componentsBufQ and
                             tbl.componentsBufQ or {1,2,3}
   -- Index in the output to store the final div(q) scalar.
   self._componentOutputDivQ = tbl.componentOutputDivQ and
                               tbl.componentOutputDivQ or 4

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.",
                        pfx, tbl._coordinate))

   self._scheme = tbl.scheme~=nil and tbl.scheme or "symmetric-node-center"
   assert(self._scheme=="symmetric-cell-center" or
          self._scheme=="symmetric-node-center",
          pfx.."scheme '"..self._scheme.."' is not supported.")

   if self._coordinate=="axisymmetric" then
      self._scheme = "symmetric-cell-center"
   end
end

local function isNaN( v ) return type( v ) == "number" and v ~= v end

function AnisotropicDiffusion:setKappaField(kappaField)
      self._kappaField = assert(self.kappaField,
                                self._pfx.."kappaMode 'field', "..
                                "Must provide 'kappaField'.")
end

-- TODO Nicer, more accurate time-step size calculation.
local suggestDt = function(kappaPara, kappaPerp, ndim, dx, cfl)
      local kappaMax = math.max(kappaPara, kappaPerp)
      local kappa__dx2_sum = 0
      for d = 1, ndim do
         kappa__dx2_sum = kappa__dx2_sum + kappaMax / ( (dx[d])^2 )
      end

      return cfl * 0.5 / kappa__dx2_sum
end

function AnisotropicDiffusion:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt

   local ndim = grid:ndim()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())
   local dx = {grid:dx(1), grid:dx(2), grid:dx(3)}

   local temp = inFld[1]
   local tempIdxr = temp:genIndexer()
   local tempPtr = temp:get(1)
   local tempPtrP = temp:get(1)
   local tempPtrM = temp:get(1)

   local emf = inFld[2]
   local emfIdxr = emf:genIndexer()
   local emfPtr = emf:get(1)

   local kappaField = inFld[3]  -- Only needed for kappamode=='field'
   local kappaFieldIdxr, kappaFieldPtr
   local useKappaField = self._kappaMode=='field'
   if (useKappaField) then
      kappaFieldIdxr = kappaField:genIndexer()
      kappaFieldPtr = kappaField:get(1)
   end

   local aux = inFld[3]  -- Only needed for kappamode=='function'
   local auxIdxr, auxPtr
   local useKappaFunction = self._kappaMode=='function'
   local hasAuxField = false
   if (useKappaFunction and type(aux)=='table') then
      auxIdxr = aux:genIndexer()
      auxPtr = aux:get(1)
      hasAuxField = true
   end

   local divQ = outFld[1]
   local divQIdxr = divQ:genIndexer()
   local divQPtr = divQ:get(1)

   local buf = outFld[2]
   local bufIdxr = buf:genIndexer()
   local bufPtr = buf:get(1)
   local bufPtrP = buf:get(1)
   local bufPtrM = buf:get(1)

   local localRange = temp:localRange()

   local status, dtSuggested = true, GKYL_MAX_DOUBLE

   local kappaPara = self._kappaPara
   local kappaPerp = self._kappaPerp
   -- Check time-step size.
   if self._kappaMode=="constant" then
      dtSuggested = suggestDt(kappaPara, kappaPerp, ndim, dx, self._cfl)
      status = dt <= dtSuggested
      if not status then return false, dtSuggested end
   end

   local cDivQ = self._componentOutputDivQ
   local cQ = self._componentsBufQ
   for icomp, c in ipairs(self._components) do

      if self._scheme=="symmetric-cell-center" then

         -- Comptue grad(T) in internal + one ghost cell centers.
         -- The resulting grad(T) components are stored in bufPtr[1,2,3].
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

               bufPtr[cQ[d]] = (tempPtrP[c] - tempPtrM[c]) * 0.5 / grid:dx(d)
            end
         end

         -- Compute q = q_para + q_perp at cell centers.
         -- q_para = -kappa_para*gradPara(T), q_perp = -kappa_perp*gradPerp(T).
         -- The resulting q components are stored in bufPtr[1,2,3] to overwrite
         -- the previously stored grad(T) components.
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

            local bDotGradT = bx*bufPtr[cQ[1]]
            if ndim>1 then bDotGradT = bDotGradT + by*bufPtr[cQ[2]] end
            if ndim>2 then bDotGradT = bDotGradT + bz*bufPtr[cQ[3]] end

            if (useKappaField) then
               kappaField:fill(kappaFieldIdxr(idx), kappaFieldPtr)
               kappaPara, kappaPerp = kappaField[1], kappaField[2]
            elseif (useKappaFunction) then
               temp:fill(tempIdxr(idx), tempPtr)
               if hasAuxField then
                  aux:fill(auxIdxr(idx), auxPtr)
               end
               kappaPara, kappaPerp = 
                  self._kappaFunction(tempPtr, emfPtr, auxPtr)
            end
            if (useKappaField or useKappaFunction) then
               dtSuggested = math.min(
                  dtSuggested,
                  suggestDt(kappaPara, kappaPerp, ndim, dx, self._cfl)
               )
               status = status and (dt <= dtSuggested)
               if not status then return false, dtSuggested end
            end

            local gradParaTx = bx * bDotGradT
            gradPerpTx = bufPtr[cQ[1]] - gradParaTx
            bufPtr[cQ[1]] = - kappaPara*gradParaTx - kappaPerp*gradPerpTx

            if ndim>1 then
               local gradParaTy = by * bDotGradT
               gradPerpTy = bufPtr[cQ[2]] - gradParaTy
               bufPtr[cQ[2]] = - kappaPara*gradParaTy - kappaPerp*gradPerpTy
            end

            if ndim>2 then
               local gradParaTz = bz * bDotGradT
               gradPerpTz = bufPtr[cQ[3]] - gradParaTz
               bufPtr[cQ[3]] = - kappaPara*gradParaTz - kappaPerp*gradPerpTz
            end
         end

         -- Compute div(q) and store it in divQPtr.
         if self._coordinate=="cartesian" then
            for idx in localRange:rowMajorIter() do
               local divq = 0
               for d = 1, ndim do
                  idx:copyInto(idxp)
                  idx:copyInto(idxm)
                  idxp[d] = idx[d]+1
                  idxm[d] = idx[d]-1

                  buf:fill(bufIdxr(idxp), bufPtrP)
                  buf:fill(bufIdxr(idxm), bufPtrM)

                  divq = divq + (bufPtrP[cQ[d]]-bufPtrM[cQ[d]]) * 0.5/grid:dx(d)
               end

               divQ:fill(divQIdxr(idx), divQPtr)
               divQPtr[cDivQ] = divq
            end
         elseif self._coordinate=="axisymmetric" then
            local xc = Lin.Vec(ndim)
            local xp = Lin.Vec(ndim)
            local xm = Lin.Vec(ndim)
            for idx in localRange:rowMajorIter() do
               local divq = 0
               for d = 1, ndim do
                  idx:copyInto(idxp)
                  idx:copyInto(idxm)
                  idxp[d] = idx[d]+1
                  idxm[d] = idx[d]-1

                  buf:fill(bufIdxr(idxp), bufPtrP)
                  buf:fill(bufIdxr(idxm), bufPtrM)

                  if d==1 then  -- R
                     grid:setIndex(idx)
                     grid:cellCenter(xc)
                     grid:setIndex(idxp)
                     grid:cellCenter(xp)
                     grid:setIndex(idxm)
                     grid:cellCenter(xm)
                     local r = xc[1]
                     local rp = xp[1]
                     local rm = xm[1]

                     divq = divq +
                            (rp*bufPtrP[cQ[d]]-rm*bufPtrM[cQ[d]]) *
                            0.5/grid:dx(d)/r
                  elseif d==2 then  -- Theta
                  elseif d==3 then  -- Z
                     divq = divq + (bufPtrP[cQ[d]]-bufPtrM[cQ[d]]) *
                            0.5/grid:dx(d)
                  end
               end

               buf:fill(bufIdxr(idx), bufPtr)
               bufPtr[5] = divq
            end
         end  -- End of different coordinate handling.

      elseif self._scheme=="symmetric-node-center" then

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
               bufPtr[cQ[d]] = (tempPtrP[c] - tempPtrM[c]) / grid:dx(d)
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
                  bufPtr[cQ[d]] = gradT * 0.5 / grid:dx(d) 
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
                  bufPtr[cQ[d]] = gradT * 0.25 / grid:dx(d) 
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

                     if (useKappaField) then
                        kappaField:fill(kappaFieldIdxr(idxp), kappaFieldPtr)
                        kappaPara, kappaPerp = kappaField[1], kappaField[2]
                     elseif (useKappaFunction) then
                        temp:fill(tempIdxr(idxp), tempPtr)
                        if hasAuxField then
                           aux:fill(auxIdxr(idx), auxPtr)
                        end
                        kappaPara, kappaPerp =
                           self._kappaFunction(tempPtr, emfPtr, auxPtr)
                     end

                     nPts = nPts + 1
                  end
               end
            end
            bx = bx / nPts
            by = by / nPts
            bz = bz / nPts

            if (useKappaField or useKappaFunction) then
               kappaPara = kappaPara / nPts
               kappaPerp = kappaPerp / nPts
               dtSuggested = math.min(
                  dtSuggested,
                  suggestDt(kappaPara, kappaPerp, ndim, dx, self._cfl)
               )
               status = status and (dt <= dtSuggested)
               if not status then return false, dtSuggested end
            end

            local bmag = math.sqrt(bx*bx + by*by + bz*bz)
            bx = bx / bmag
            by = by / bmag
            bz = bz / bmag
            assert(bmag>0, "Zero B field detected!")

            -- Compute gradParaT and gradPerpT at cell corners.
            buf:fill(bufIdxr(idx), bufPtr)
            local bDotGradT = bx*bufPtr[cQ[1]]+by*bufPtr[cQ[2]]+bz*bufPtr[cQ[3]]

            local gradParaTx = bx * bDotGradT
            local gradParaTy = by * bDotGradT
            local gradParaTz = bz * bDotGradT

            gradPerpTx = bufPtr[cQ[1]] - gradParaTx
            gradPerpTy = bufPtr[cQ[2]] - gradParaTy
            gradPerpTz = bufPtr[cQ[3]] - gradParaTz

            buf:fill(bufIdxr(idx), bufPtr)

            bufPtr[cQ[1]] = - kappaPara*gradParaTx - kappaPerp*gradPerpTx
            bufPtr[cQ[2]] = - kappaPara*gradParaTy - kappaPerp*gradPerpTy
            bufPtr[cQ[3]] = - kappaPara*gradParaTz - kappaPerp*gradPerpTz
         end  -- Node-center q computation ends.

         -- Compute div(q) at cell centers.
         for idx in localRange:rowMajorIter() do
            divQ:fill(divQIdxr(idx), divQPtr)
            local divq = 0

            if ndim==1 then
               local d = 1
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]
               buf:fill(bufIdxr(idxp), bufPtrP)
               buf:fill(bufIdxr(idxm), bufPtrM)
               divq = divq + (bufPtrP[cQ[d]]-bufPtrM[cQ[d]]) / grid:dx(d)
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
                     divq = divq+(bufPtrP[cQ[d]]-bufPtrM[cQ[d]])*0.5/grid:dx(d)
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
                        divq = divq + (bufPtrP[cQ[d]]-bufPtrM[cQ[d]]) *
                               0.25/grid:dx(d)
                     end
                  end
               end -- Loop over gradient directions.
            end -- Loop over ndim==1,2,3 ends.
            divQPtr[cDivQ] = divq 
         end -- div(q) computation ends.

      end  -- scheme handling in divq calculation ends.

      -- Add div(q) to energy.
      for idx in localRange:rowMajorIter() do
         temp:fill(tempIdxr(idx), tempPtr)
         divQ:fill(divQIdxr(idx), divQPtr)
         tempPtr[c] = tempPtr[c] - dt * divQPtr[cDivQ]
      end

   end -- Ends of loop over components of which diffusion if computed.

   return status, dtSuggested
end

function AnisotropicDiffusion:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return AnisotropicDiffusion
