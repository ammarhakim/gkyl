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

   self._mass = assert(tbl.mass, pfx .. "Must provide 'mass'.")

   self._charge = assert(tbl.charge, pfx .. "Must provide 'charge'.")
   assert(#self._mass==self._nFluids and #self._charge==self._nFluids,
          pfx .. " Lengths of mass of charge must match nFluids.")

   self._epsilon0 = assert(tbl.epsilon0, pfx .. "Must provide 'epsilon0'.")

   self._logA = tbl.coulombLogarithm and tbl.coulombLogarithm or 10

   -- Calculate tau_e with plasma parameters vs. using a preset value
   self._calcTau = tbl.calcTau ~= nil and tbl.calcTau or false
   self._tau = tbl.tau
   assert(not (type(self._tau)=='table' and self.calcTau),
          pfx ..  "Cannot specify 'tau' and 'calcTau' simultaneously.")
   assert(type(self._tau)=='table' or self.calcTau,
          pfx ..  "Must specify one of 'tau' and 'calcTau'.")

   self._hasHeating = tbl.hasHeating ~= nil and tbl.hasHeating or false

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.", pfx, tbl._coordinate))

   assert(self._gasGamma==5./3., pfx .. " gasGamma must be 5/3.")
end

local pressure = function (q, gasGamma, mass)
	 return (gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
end

local temperature = function (q, gasGamma, mass)
	 return pressure(q, gasGamma, mass) * mass / q[1]
end

local calcTau = function(fluid, mass, charge)
   local tau = 1 -- TODO
   return tau
end

local calcEta = function(fluidPtr, emfPtr, mass, charge, gasGamma)
   local tau = calcTau(fluidPtr, mass, charge)
   local pr = pressure(fluidPtr, gasGamma, mass)
   local bx = emfPtr[4]
   local by = emfPtr[5]
   local bz = emfPtr[6]
   local bmag = math.sqrt(bx*bx + by*by + bz*bz)
   local Omega = math.abs(charge*bmag/mass)
   local eta = pr * tau / (Omega*tau^2)
   return eta
end

function BraginskiiViscosityDiffusion:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local gasGamma = self._gasGamma
   local epsilon0 = self._epsilon0
   local nFluids = self._nFluids
   local logA = self._logA

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   -- Indices.
   local ndim = grid:ndim()
   local idxp = Lin.IntVec(grid:ndim())
   local idxm = Lin.IntVec(grid:ndim())

   local emf = outFld[nFluids+1]
   local emfIdxr = emf:genIndexer()
   local emfPtr = emf:get(1)
   local emfPtrP = emf:get(1)
   local emfPtrM = emf:get(1)

   -- Comptue grad_para(T) ain internal cells.
   for s = 1, nFluids do
      local fluid = outFld[s]
      local fluidIdxr = fluid:genIndexer()
      local fluidPtr = fluid:get(1)
      local fluidPtrP = fluid:get(1)
      local fluidPtrM = fluid:get(1)

      local fluidBuf = inFld[s]
      local fluidBufIdxr = fluidBuf:genIndexer()
      local fluidBufPtr = fluidBuf:get(1)
      local fluidBufPtrP = fluidBuf:get(1)
      local fluidBufPtrM = fluidBuf:get(1)

      local mass = self._mass[s]
      local charge = self._charge[s]

      local localRange = fluid:localRange()

      -- Compute rhs.
      if self._coordinate=="cartesian" then
         for idx in localRange:rowMajorIter() do
            local qVis = 0
            for d = 1, ndim do
               idx:copyInto(idxp)
               idx:copyInto(idxm)
               idxp[d] = idx[d]+1
               idxm[d] = idx[d]-1
              
               fluid:fill(fluidIdxr(idx), fluidPtr)
               fluid:fill(fluidIdxr(idxp), fluidPtrP)
               fluid:fill(fluidIdxr(idxm), fluidPtrM)
               emf:fill(emfIdxr(idxp), emfPtrP)
               emf:fill(emfIdxr(idxm), emfPtrM)

               local eta = calcEta(fluidPtr, emfPtr, mass, charge, gasGamma)
               local etaP = calcEta(fluidPtrP, emfPtrP, mass, charge, gasGamma)
               local etaM = calcEta(fluidPtrM, emfPtrM, mass, charge, gasGamma)

               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               -- Compute momentum diffusion as grad(eta * grad(V)).
               local dx = grid:dx(d)
               for c=2,4 do
                  fluidBufPtr[c] = (etaP*(fluidPtrP[c]-fluidPtr[c]) -
                                    etaM*(fluidPtr[c]-fluidPtrM[c])) / (dx*dx)
               end

               -- Compute viscous heating eta*grad(V):grad(V) as
               -- eta*|grad(V)|^2.
               if self._hasHeating then
                  local eta = calcEta(fluidPtr, emfPtr, mass, charge, gasGamma, s)
                  for c=2,4 do
                     qVis = qVis +
                            eta * ( (fluidPtrP[c]-fluidPtrM[c])/(2*dx) )^2
                  end
               end
            end
            fluidBufPtr[5] = qVis
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
               fluid:fill(fluidIdxr(idx), fluidPtr)
               fluid:fill(fluidIdxr(idxp), fluidPtrP)
               fluid:fill(fluidIdxr(idxm), fluidPtrM)
               emf:fill(emfIdxr(idxp), emfPtrP)
               emf:fill(emfIdxr(idxm), emfPtrM)
               local eta = calcEta(fluidPtr, emfPtr, mass, charge, gasGamma)
               local etaP = calcEta(fluidPtrP, emfPtrP, mass, charge, gasGamma)
               local etaM = calcEta(fluidPtrM, emfPtrM, mass, charge, gasGamma)

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

               fluidBufPtr[2] = (etaPH*rpH*(fluidPtrP[2]-fluidPtr [2]) -
                                 etaMH*rmH*(fluidPtr [2]-fluidPtrM[2])) / (dr*dr*r) -
                                eta*fluidPtr[2]/(r*r)
               fluidBufPtr[3] = (etaPH*(rpH^3)*(fluidPtrP[3]/rp-fluidPtr [3]/r) -
                                 etaMH*(rmH^3)*(fluidPtr [3]/r -fluidPtrM[3]/rm)) / (dr*dr*r)
               fluidBufPtr[4] = (etaPH*rpH*(fluidPtrP[4]-fluidPtr [4]) -
                                 etaMH*rmH*(fluidPtr [4]-fluidPtrM[4])) / (dr*dr*r)

               if self._hasHeating then
                  local eta = calcEta(fluidPtr, emfPtr, mass, charge, gasGamma)
                  qVis = qVis +
                         eta * (fluidPtrP[3]/rp-fluidPtrM[3]/rm) * 0.5 / dr
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
               fluid:fill(fluidIdxr(idx), fluidPtr)
               fluid:fill(fluidIdxr(idxp), fluidPtrP)
               fluid:fill(fluidIdxr(idxm), fluidPtrM)
               emf:fill(emfIdxr(idxp), emfPtrP)
               emf:fill(emfIdxr(idxm), emfPtrM)
               local eta = calcEta(fluidPtr, emfPtr, mass, charge, gasGamma)
               local etaP = calcEta(fluidPtrP, emfPtrP, mass, charge, gasGamma)
               local etaM = calcEta(fluidPtrM, emfPtrM, mass, charge, gasGamma)

               local etaPH = 0.5 * (eta+etaP)
               local etaMH = 0.5 * (eta+etaM)

               for c=2,4 do
                  fluidBufPtr[c] = fluidBufPtr[c] +
                                   (etaPH*(fluidPtrP[c]-fluidPtr [c]) -
                                    etaMH*(fluidPtr[c] -fluidPtrM[c])) / (dz*dz)
               end

               if self._hasHeating then
                  grid:setIndex(idx)
                  grid:cellCenter(xc)
                  qVis = qVis +
                         eta * (fluidPtrP[3]-fluidPtrM[3]) / r * 0.5 / dz
               end
            end

            fluidBufPtr[5] = qVis
         end
      else
         assert(false)
      end
   end

   -- Apply rhs.
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

         if self._hasHeating then
            local keOld = 0.5*(fluidPtr[2]^2+fluidPtr[3]^2+fluidPtr[4]^2)/fluidPtr[1]
            for c=2,4 do
               fluidPtr[c] = fluidPtr[c] + fluidBufPtr[c]
            end
            local keNew = 0.5*(fluidPtr[2]^2+fluidPtr[3]^2+fluidPtr[4]^2)/fluidPtr[1]
            fluidPtr[5] = fluidPtr[5]+keNew-keOld+fluidBufPtr[5]
         else
            for c=2,4 do
               fluidPtr[c] = fluidPtr[c] + fluidBufPtr[c]
            end
         end
      end
   end

   return status, dtSuggested
end

function BraginskiiViscosityDiffusion:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiViscosityDiffusion
