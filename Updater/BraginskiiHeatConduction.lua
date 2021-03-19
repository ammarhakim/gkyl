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

   local heatFlux = inFld[nFluids+1]  -- This is really emfBuf.
   local heatFluxIdxr = heatFlux:genIndexer()
   local heatFluxPtr = heatFlux:get(1)
   local heatFluxPtrP = heatFlux:get(1)
   local heatFluxPtrM = heatFlux:get(1)

   for s = 1, nFluids do
      local fluid = outFld[s]
      local fluidIdxr = fluid:genIndexer()
      local fluidPtr = fluid:get(1)

      local fluidBuf = inFld[s]
      local fluidBufIdxr = fluidBuf:genIndexer()
      local fluidBufPtr = fluidBuf:get(1)
      local fluidBufPtrP = fluidBuf:get(1)
      local fluidBufPtrM = fluidBuf:get(1)

      local mass = self._mass[s]
      local charge = self._charge[s]

      local localExtRange = fluid:localExtRange()
      local localRange = fluid:localRange()
      local localExt1Range = localRange:extend(1, 1)

      -- Compute temperature in internal and ghost cell centers.
      for idx in localExtRange:rowMajorIter() do
         fluid:fill(fluidIdxr(idx), fluidPtr)
         fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)

         fluidBufPtr[5] = temperature(fluidPtr, gasGamma, mass)
      end

      -- Comptue grad(T) in internal + one ghost cell centers.
      for idx in localExt1Range:rowMajorIter() do
         for d = 1, ndim do
            idx:copyInto(idxp)
            idx:copyInto(idxm)
            idxp[d] = idx[d]+1
            idxm[d] = idx[d]-1

            fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)
            fluidBuf:fill(fluidBufIdxr(idxp), fluidBufPtrP)
            fluidBuf:fill(fluidBufIdxr(idxm), fluidBufPtrM)
            fluidBufPtr[d+1] = (fluidBufPtrP[5] - fluidBufPtrM[5]) *
                               0.5 / grid:dx(d)
         end
      end

      -- Compute q = q_para + q_perp at cell centers.
      -- q_para = kappa_para*gradPara(T), q_perp = kappa_perp*gradPerp(T).
      -- For the electron fluid in a two-fluid plasma, also add -0.71*pe*dVpara.
      for idx in localExt1Range:rowMajorIter() do
         emf:fill(emfIdxr(idx), emfPtr)
         fluid:fill(fluidIdxr(idx), fluidPtr)
         fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)
         heatFlux:fill(heatFluxIdxr(idx), heatFluxPtr)

         local bx = emfPtr[4]
         local by = emfPtr[5]
         local bz = emfPtr[6]
         local bmag = math.sqrt(bx*bx + by*by + bz*bz)
         bx = bx / bmag
         by = by / bmag
         bz = bz / bmag
         assert(bmag>0, "Zero B field detected!")

         local bDotGradT = bx*fluidBufPtr[2]+by*fluidBufPtr[3]+bz*fluidBufPtr[4]

         local gradParaTx = bx * bDotGradT
         local gradParaTy = by * bDotGradT
         local gradParaTz = bz * bDotGradT

         gradPerpTx = fluidBufPtr[2] - gradParaTx
         gradPerpTy = fluidBufPtr[3] - gradParaTy
         gradPerpTz = fluidBufPtr[4] - gradParaTz

         local n = fluidPtr[1] / mass
         local T = fluidBufPtr[5]

         -- TODO: Calculate tau.
         local tau = self._tau[s]
         local Omega = math.abs(charge*bmag)/mass
         -- TODO Provide constant coefficients.
         local kappaPara = -n*T*tau/mass
         local kappaPerp = -n*T*tau/mass/(Omega*tau)^2

         heatFluxPtr[1] = kappaPara*gradParaTx + kappaPerp*gradPerpTx
         heatFluxPtr[2] = kappaPara*gradParaTy + kappaPerp*gradPerpTy
         heatFluxPtr[3] = kappaPara*gradParaTz + kappaPerp*gradPerpTz

         -- FIXME: Nicer handling of terms that involve other species.
         if nFluids==2 and charge<0 then
            local elcPtr = fluidPtr
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
            heatFluxPtr[1] = heatFluxPtr[1] + 0.71*pr*dVparx
            heatFluxPtr[2] = heatFluxPtr[2] + 0.71*pr*dVpary
            heatFluxPtr[3] = heatFluxPtr[3] + 0.71*pr*dVparz
         end
      end

      -- Compute div(q) at cell-centers using cell-center q values.
      if self._coordinate=="cartesian" then
         for idx in localRange:rowMajorIter() do
            fluid:fill(fluidIdxr(idx), fluidPtr)
            fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)
            heatFlux:fill(heatFluxIdxr(idx), heatFluxPtr)

            local divq = 0
            for d = 1, ndim do
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1

               heatFlux:fill(heatFluxIdxr(idxp), heatFluxPtrP)
               heatFlux:fill(heatFluxIdxr(idxm), heatFluxPtrM)

               divq = divq + (heatFluxPtrP[d] - heatFluxPtrM[d]) *
                      0.5 / grid:dx(d)
            end

           fluidBufPtr[1] = divq
         end
      elseif self._coordinate=="axisymmetric" then
         local xc = Lin.Vec(ndim)
         local xp = Lin.Vec(ndim)
         local xm = Lin.Vec(ndim)
         for idx in localRange:rowMajorIter() do
            fluid:fill(fluidIdxr(idx), fluidPtr)
            fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)
            heatFlux:fill(heatFluxIdxr(idx), heatFluxPtr)

            local divq = 0
            for _,d in ipairs({1,3}) do
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1

               heatFlux:fill(heatFluxIdxr(idxp), heatFluxPtrP)
               heatFlux:fill(heatFluxIdxr(idxm), heatFluxPtrM)

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
                         (rp*heatFluxPtrP[d]-rm*heatFluxPtrM[d]) *
                         0.5/grid:dx(d)/r
               elseif d==2 then  -- Theta
               elseif d==3 then  -- Z
                  divq = divq +
                         (heatFluxPtrP[d]-heatFluxPtrM[d]) * 0.5/grid:dx(d)
               else
                  assert(false)
               end
            end

           fluidBufPtr[1] = divq
         end
      end
   end

   -- Add div(q) to energy.
   for s = 1, nFluids do
      local fluid = outFld[s]
      local fluidIdxr = fluid:genIndexer()
      local fluidPtr = fluid:get(1)

      local fluidBuf = inFld[s]
      local fluidBufIdxr = fluidBuf:genIndexer()
      local fluidBufPtr = fluidBuf:get(1)

      local localRange = fluid:localRange()

      for idx in localRange:rowMajorIter() do
         fluid:fill(fluidIdxr(idx), fluidPtr)
         fluidBuf:fill(fluidBufIdxr(idx), fluidBufPtr)
         fluidPtr[5] = fluidPtr[5] - fluidBufPtr[5]
      end
   end

   return status, dtSuggested
end

function BraginskiiHeatConduction:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiHeatConduction
