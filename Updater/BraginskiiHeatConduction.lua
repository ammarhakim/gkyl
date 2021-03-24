-- Gkyl ------------------------------------------------------------------------
--
-- Updater to apply Braginskii-like heat conduction.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local BraginskiiHeatConduction = Proto(UpdaterBase)

function BraginskiiHeatConduction:init(tbl)
   BraginskiiHeatConduction.super.init(self, tbl)

   local pfx = "Updater.BraginskiiHeatConduction: "

   self._onGrid = assert(tbl.onGrid, pfx.."Must provide 'onGrid'.")

   self._gasGamma = assert(tbl.gasGamma, pfx .. "Must provide 'gasGamma'.")

   self._nFluids = assert(tbl.numFluids, pfx .. "Must provide 'numFluids'.")

   self._mass = assert(tbl.mass, pfx .. "Must provide 'mass'.")

   self._charge = assert(tbl.charge, pfx .. "Must provide 'charge'.")
   assert(#self._mass==self._nFluids and #self._charge==self._nFluids,
          pfx .. "Lengths of mass of charge must match nFluids.")

   self._epsilon0 = assert(tbl.epsilon0, pfx .. "Must provide 'epsilon0'.")

   self._logA = tbl.coulombLogarithm and tbl.coulombLogarithm or 10

   -- Calculate tau_e with plasma parameters vs. using a preset value
   self._calcTau = tbl.calcTau ~= nil and tbl.calcTau or false
   self._tau = tbl.tau
   assert(not (type(self._tau)=='table' and self.calcTau),
          pfx ..  "Cannot specify 'tau' and 'calcTau' simultaneously.")
   assert(type(self._tau)=='table' or self.calcTau,
          pfx ..  "Must specify one of 'tau' and 'calcTau'.")

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.",
                        pfx, tbl._coordinate))

   self._timeStepper = tbl.timeStepper~=nil and tbl.timeStepper or
                       "two-step-node-center"
   assert(self._timeStepper=="two-step-cell-center" or
          self._timeStepper=="two-step-node-center",
          pfx.."timeStepper '"..self._timeStepper.."' is not supported.")

   assert(self._gasGamma==5./3., pfx .. "gasGamma must be 5/3.")
end

local temperature = function (q, gasGamma, mass)
	 local pr = (gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
    return pr * mass / q[1]
end

function BraginskiiHeatConduction:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local gasGamma = self._gasGamma
   local epsilon0 = self._epsilon0
   local nFluids = self._nFluids
   local logA = self._logA

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   local ndim = grid:ndim()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())

   local emf = outFld[nFluids+1]
   local emfIdxr = emf:genIndexer()
   local emfPtr = emf:get(1)

   for s = 1, nFluids do

      local fld = outFld[s]
      local fldIdxr = fld:genIndexer()
      local fldPtr = fld:get(1)

      local buf = inFld[s]
      local bufIdxr = buf:genIndexer()
      local bufPtr = buf:get(1)
      local bufPtrP = buf:get(1)
      local bufPtrM = buf:get(1)

      local mass = self._mass[s]
      local charge = self._charge[s]

      local localRange = fld:localRange()

      -- Compute temperature in internal and ghost cell centers.
      local localExtRange = fld:localExtRange()
      for idx in localExtRange:rowMajorIter() do
         fld:fill(fldIdxr(idx), fldPtr)
         buf:fill(bufIdxr(idx), bufPtr)
         bufPtr[4] = temperature(fldPtr, gasGamma, mass)
      end

      -- Different algorithms to compute q and div(q).
      if self._timeStepper=='two-step-cell-center' then

         -- Comptue grad(T) in internal + one ghost cell centers.
         local localExt1Range = localRange:extend(1, 1)
         for idx in localExt1Range:rowMajorIter() do
            for d = 1, ndim do
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1

               buf:fill(bufIdxr(idx), bufPtr)
               buf:fill(bufIdxr(idxp), bufPtrP)
               buf:fill(bufIdxr(idxm), bufPtrM)
               bufPtr[d] = (bufPtrP[4] - bufPtrM[4]) * 0.5 / grid:dx(d)
            end
         end

         -- Compute q = q_para + q_perp at cell centers.
         -- q_para = kappa_para*gradPara(T), q_perp = kappa_perp*gradPerp(T).
         -- For the electron fld in a two-fld plasma, also add
         -- -0.71*pe*dVpara.
         for idx in localExt1Range:rowMajorIter() do
            emf:fill(emfIdxr(idx), emfPtr)
            fld:fill(fldIdxr(idx), fldPtr)
            buf:fill(bufIdxr(idx), bufPtr)

            local bx = emfPtr[4]
            local by = emfPtr[5]
            local bz = emfPtr[6]
            local bmag = math.sqrt(bx*bx + by*by + bz*bz)
            bx = bx / bmag
            by = by / bmag
            bz = bz / bmag
            assert(bmag>0, "Zero B field detected!")

            local n = fldPtr[1] / mass
            local T = bufPtr[4]
            -- TODO: Calculate tau.
            local tau = self._tau[s]
            local Omega = math.abs(charge*bmag)/mass
            -- TODO Provide constant coefficients.
            local kappaPara = -n*T*tau/mass
            local kappaPerp = -n*T*tau/mass/(Omega*tau)^2

            local bDotGradT = bx*bufPtr[1]
            if ndim>1 then bDotGradT = bDotGradT + by*bufPtr[2] end 
            if ndim>2 then bDotGradT = bDotGradT + bz*bufPtr[3] end 

            local gradParaTx = bx * bDotGradT
            local gradPerpTx = bufPtr[1] - gradParaTx
            bufPtr[1] = kappaPara*gradParaTx + kappaPerp*gradPerpTx

            if ndim>1 then
               local gradParaTy = by * bDotGradT
               local gradPerpTy = bufPtr[2] - gradParaTy
               bufPtr[2] = kappaPara*gradParaTy + kappaPerp*gradPerpTy
            end

            if ndim>2 then
               local gradParaTz = bz * bDotGradT
               local gradPerpTz = bufPtr[3] - gradParaTz
               bufPtr[3] = kappaPara*gradParaTz + kappaPerp*gradPerpTz
            end

            -- FIXME: Nicer handling of terms that involve other species.
            if nFluids==2 and charge<0 then
               local elcPtr = fldPtr
               local ion = outFld[2]
               local ionIdxr = ion:genIndexer()
               local ionPtr = ion:get(1)
               ion:fill(ionIdxr(idx), ionPtr)

               local dVx = ionPtr[2]/ionPtr[1] - elcPtr[2]/elcPtr[1]
               local dVy = ionPtr[3]/ionPtr[1] - elcPtr[3]/elcPtr[1]
               local dVz = ionPtr[4]/ionPtr[1] - elcPtr[4]/elcPtr[1]
               local bDotDV = bx * dVx + by * dVy + bz * dVz
               local dVparx = bx * bDotDV
               local dVpary = by * bDotDV
               local dVparz = bz * bDotDV

               local pr = n * T
               bufPtr[1] = bufPtr[1] + 0.71*pr*dVparx
               bufPtr[2] = bufPtr[2] + 0.71*pr*dVpary
               bufPtr[3] = bufPtr[3] + 0.71*pr*dVparz
            end
         end

         -- Compute div(q) at cell-centers using cell-center q values.
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

                  divq = divq + (bufPtrP[d]-bufPtrM[d])*0.5/grid:dx(d)
               end

               buf:fill(bufIdxr(idx), bufPtr)
               bufPtr[5] = divq
            end
         elseif self._coordinate=="axisymmetric" then
            local xc = Lin.Vec(ndim)
            local xp = Lin.Vec(ndim)
            local xm = Lin.Vec(ndim)
            for idx in localRange:rowMajorIter() do
               local divq = 0
               for _,d in ipairs({1,3}) do
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

                     divq = divq+(rp*bufPtrP[d]-rm*bufPtrM[d])*0.5/grid:dx(d)/r
                  elseif d==2 then  -- Theta
                  elseif d==3 then  -- Z
                     divq = divq + (bufPtrP[d]-bufPtrM[d])*0.5/grid:dx(d)
                  else
                     assert(false)
                  end
               end

               buf:fill(bufIdxr(idx), bufPtr)
               bufPtr[5] = divq
            end
         end

      elseif self._timeStepper=='two-step-node-center' then

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
               buf:fill(bufIdxr(idxp), bufPtrP)
               buf:fill(bufIdxr(idxm), bufPtrM)
               bufPtr[d] = (bufPtrP[4] - bufPtrM[4]) / grid:dx(d)
            elseif ndim==2 then
               local subIterDirs = {{2}, {1}}
               for d=1,ndim do  -- Gradient direction.
                  local gradT = 0

                  idx:copyInto(idxp)
                  idx:copyInto(idxm)
                  idxp[d] = idx[d]
                  idxm[d] = idx[d]-1

                  -- Add contributions from cell-center values of cells sharing
                  -- this corner.
                  local i1 = subIterDirs[d][1]
                  for _,commonShift1 in ipairs({-1, 0}) do
                     idxp[i1] = idx[i1] + commonShift1
                     idxm[i1] = idx[i1] + commonShift1
                     buf:fill(bufIdxr(idxp), bufPtrP)
                     buf:fill(bufIdxr(idxm), bufPtrM)
                     gradT = gradT + (bufPtrP[4] - bufPtrM[4])
                  end
                  bufPtr[d] = gradT * 0.5 / grid:dx(d) 
               end
            elseif ndim==3 then
               local subIterDirs = {{2,3}, {1,3}, {1,2}}
               for d=1,ndim do  -- Gradient direction.
                  local gradT = 0

                  idx:copyInto(idxp)
                  idx:copyInto(idxm)
                  idxp[d] = idx[d]
                  idxm[d] = idx[d]-1

                  -- Add contributions from cell-center values of cells sharing
                  -- this corner.
                  local i1 = subIterDirs[d][1]
                  local i2 = subIterDirs[d][2]
                  for _,commonShift1 in ipairs({-1, 0}) do
                     for _,commonShift2 in ipairs({-1, 0}) do
                        idxp[i1] = idx[i1] + commonShift1
                        idxm[i1] = idx[i1] + commonShift1
                        idxp[i2] = idx[i2] + commonShift2
                        idxm[i2] = idx[i2] + commonShift2
                        buf:fill(bufIdxr(idxp), bufPtrP)
                        buf:fill(bufIdxr(idxm), bufPtrM)
                        gradT = gradT + (bufPtrP[4] - bufPtrM[4])
                     end
                  end
                  bufPtr[d] = gradT * 0.25 / grid:dx(d) 
               end
            end -- Loop over ndim==1,2,3 ends.
         end -- Node-center grad(T) computation ends.

         -- Compute q on nodes (cell corners).
         -- Note that fld values and buf[4] T values are at cell centers,
         -- buf grad(T) values and heat flux values are at cell corners.
         for idx in localNodeRange:rowMajorIter() do
            -- Compute B field, n and T (for kappa) at cell corners.
            -- The value at the i-th node is defined as an average of values at
            -- centers of all cells sharing this node.
            local bx, by, bz = 0, 0, 0
            local n = 0
            local T = 0

            local nPts = 1
            local xshifts = {-1, 0}
            local yshifts = ndim>1 and {-1, 0} or {}
            local zshifts = ndim>2 and {-1, 0} or {}
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

                     fld:fill(fldIdxr(idxp), fldPtr)
                     n = n + fldPtr[1] / mass

                     -- T values are stored at cell-centers.
                     buf:fill(bufIdxr(idxp), bufPtr)
                     T = T + bufPtr[4]

                     nPts = nPts + 1
                  end
               end
            end
            bx = bx / nPts
            by = by / nPts
            bz = bz / nPts
            n = n / nPts
            T = T / nPts

            local bmag = math.sqrt(bx*bx + by*by + bz*bz)
            bx = bx / bmag
            by = by / bmag
            bz = bz / bmag
            assert(bmag>0, "Zero B field detected!")

            -- Compute kappaPara and kappaPerp at cell corners.
            -- TODO: Calculate tau.
            local tau = self._tau[s]
            local Omega = math.abs(charge*bmag)/mass
            -- TODO Provide constant coefficients.
            local kappaPara = -n*T*tau/mass
            local kappaPerp = -n*T*tau/mass/(Omega*tau)^2

            -- Compute gradParaT and gradPerpT at cell corners.
            buf:fill(bufIdxr(idx), bufPtr)
            local bDotGradT = bx*bufPtr[1]
            if ndim>1 then bDotGradT = bDotGradT + by*bufPtr[2] end 
            if ndim>2 then bDotGradT = bDotGradT + bz*bufPtr[3] end 

            -- Compute heat flux at cell corners.
            local gradParaTx = bx * bDotGradT
            gradPerpTx = bufPtr[1] - gradParaTx
            bufPtr[1] = kappaPara*gradParaTx + kappaPerp*gradPerpTx

            if ndim>1 then
               local gradParaTy = by * bDotGradT
               gradPerpTy = bufPtr[2] - gradParaTy
               bufPtr[2] = kappaPara*gradParaTy + kappaPerp*gradPerpTy
            end

            if ndim>2 then
               local gradParaTz = bz * bDotGradT
               gradPerpTz = bufPtr[3] - gradParaTz
               bufPtr[3] = kappaPara*gradParaTz + kappaPerp*gradPerpTz
            end
         end  -- Node-center q computation ends.

         -- Compute div(q) at cell centers.
         for idx in localRange:rowMajorIter() do

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
                  local i1 = subIterDirs[d][1]
                  for _,commonShift1 in ipairs({0, 1}) do
                     idxp[i1] = idx[i1] + commonShift1
                     idxm[i1] = idx[i1] + commonShift1
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
                  local i1 = subIterDirs[d][1]
                  local i2 = subIterDirs[d][2]
                  for _,commonShift1 in ipairs({0, 1}) do
                     for _,commonShift2 in ipairs({0, 1}) do
                        idxp[i1] = idx[i1] + commonShift1
                        idxm[i1] = idx[i1] + commonShift1
                        idxp[i2] = idx[i2] + commonShift2
                        idxm[i2] = idx[i2] + commonShift2
                        buf:fill(bufIdxr(idxp), bufPtrP)
                        buf:fill(bufIdxr(idxm), bufPtrM)
                        divq = divq + (bufPtrP[d]-bufPtrM[d])*0.25/grid:dx(d)
                     end
                  end
               end
            end -- Different handling of 1d/2d/3d n div(q) computation ends.

            buf:fill(bufIdxr(idx), bufPtr)
            bufPtr[5] = divq 
         end -- div(q) computation ends.

      end  -- timeStepper handling in divq calculation ends.

   end  -- Loop over species ends.

   -- Add div(q) to energy.
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
         fldPtr[5] = fldPtr[5] - bufPtr[5]
      end
   end

   return status, dtSuggested
end

function BraginskiiHeatConduction:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiHeatConduction
