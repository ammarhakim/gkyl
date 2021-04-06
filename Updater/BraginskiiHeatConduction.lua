-- Gkyl ------------------------------------------------------------------------
--
-- Updater to apply Braginskii-like heat conductivity to the fluid energy
-- equation.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local AnisotropicDiffusion = require "Updater.AnisotropicDiffusion"
local xsys = require "xsys"

local BraginskiiHeatConduction = Proto(UpdaterBase)

function BraginskiiHeatConduction:init(tbl)
   BraginskiiHeatConduction.super.init(self, tbl)
   local pfx = "Updater.BraginskiiHeatConduction: "

   self._onGrid = assert(tbl.onGrid, pfx.."Must provide 'onGrid'.")
   self._cfl = tbl.cfl~=nil and tbl.cfl or 0.9

   self._gasGamma = assert(tbl.gasGamma, pfx .. "Must provide 'gasGamma'.")

   self._nFluids = assert(tbl.numFluids, pfx .. "Must provide 'numFluids'.")
   self._mass = assert(tbl.mass, pfx .. "Must provide 'mass'.")
   self._charge = assert(tbl.charge, pfx .. "Must provide 'charge'.")
   assert(#self._mass==self._nFluids and #self._charge==self._nFluids,
          pfx.."Number of values in 'mass' or 'charge'  must match 'nFluids'.")

   self._kappaMode = tbl.kappaMode and tbl.kappaMode or "constant"
   if self._kappaMode=="constant" then
      local pfxx = pfx.."For kappaMode=='constant', "
      self._kappaPara = assert(tbl.kappaPara, pfxx.."Must provide 'kappaPara'.")
      self._kappaPerp = assert(tbl.kappaPerp, pfxx.."Must provide 'kappaPerp'.")

      assert(#self._kappaPara==self._nFluids,
             pfxx .. "'kappaPara' must be a list with 'nFluids' entries.")
      assert(#self._kappaPerp==self._nFluids,
             pfxx .. "'kappaPerp' must be a list with 'nFluids' entries.")
      for s=1,self._nFluids do
         assert(self._kappaPara[s]>=0, pfxx.."kappaPara values must be >=0.")
         assert(self._kappaPerp[s]>=0, pfxx.."kappaPerp values must be >=0.")
      end
   elseif self._kappaMode=="from-tau" then
      local pfxx = pfx.."For kappaMode=='from-tau', "
      self._tau = tbl.tau
      assert(type(self._tau)=="table" and #self._tau==self._nFluids,
             pfx .. "For kappaMode=='from-tau', 'tau' must be a list with "..
             "'nFluids' entries.")
   elseif self._kappaMode=="from-scratch" then
      assert(false, pfx.."kappaMode 'from-scratch' is not implemented.")
   else
      assert(false, pfx.."kappaMode "..self._kappaMode.." is not recognized.")
   end

   self._coordinate = tbl.coordinate ~= nil and tbl.coordinate or "cartesian"
   assert(self._coordinate=="cartesian" or
          self._coordinate=="axisymmetric",
          string.format("%s'coordinate' %s not recognized.",
                        pfx, tbl._coordinate))

   self._scheme = tbl.scheme~=nil and tbl.scheme or "symmetric-cell-center"
   assert(self._scheme=="symmetric-cell-center" or
          self._scheme=="symmetric-node-center",
          pfx.."scheme '"..self._scheme.."' is not supported.")

   if self._coordinate=="axisymmetric" then
      self._scheme = "symmetric-cell-center"
   end

   local anisotropicDiffusionKappaMode
   if self._kappaMode=="constant" then
      anisotropicDiffusionKappaMode = "constant"
   elseif self.._kappaMode=="from-tau" then
      anisotropicDiffusionKappaMode = "function"
   elseif self.._kappaMode=="from-scratch" then
      anisotropicDiffusionKappaMode = "function"
   end
   self._anisotropicDiffusion = AnisotropicDiffusion {
      onGrid = self._onGrid,
      cfl = self._cfl,
      scheme = self._scheme,
      coordinate = self._coordinate,
      components = {4},  -- component in buffer to store the temperature, T
      componentsBufQ = {1,2,3},  -- components in buffer to store q vector
      componentOutputDivQ = 5,  -- component in buffer to store div(q)
      updateTemperature = false,

      kappaMode = anisotropicDiffusionKappaMode,
      kappaPara = 1,
      kappaPerp = 1,
   }

   -- Add additional Braginskii terms. FIXME Not implemented yet.
   self._addAdditionalQ = xsys.pickBool(tbl.addAdditionalQ, false)
   assert(not self._addAdditionalQ,
          pfx.."addAdditionalQ option is not implemented.")
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
   local nFluids = self._nFluids

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   local ndim = grid:ndim()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())

   local emf = outFld[nFluids+1]
   local emfIdxr = emf:genIndexer()
   local emfPtr = emf:get(1)

   local staticEmf = inFld[nFluids+2]

   local localRange = emf:localRange()
   local localExtRange = emf:localExtRange()
   local localExt1Range = localRange:extend(1, 1)

   self._anisotropicDiffusion:setDtAndCflRate(myDt)

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

      -- Step 1: Compute temperature in internal and ghost cell centers.
      for idx in localExtRange:rowMajorIter() do
         fld:fill(fldIdxr(idx), fldPtr)
         buf:fill(bufIdxr(idx), bufPtr)
         bufPtr[4] = temperature(fldPtr, gasGamma, mass)
      end

      -- Step 2: Compute q.
      if self._kappaMode=="constant" then
         self._anisotropicDiffusion:setKappa(self._kappaPara[s],
                                             self._kappaPerp[s])
      elseif self._kappaMode=="from-tau" then
         local tau = self._tau[s]
         local kappaFunction = function(
               bmag, bufPtr, emfPtr, fldPtr, staticEmfPtr)
            local n = fldPtr[1] / mass
            local T = bufPtr[4]
            local Omega = math.abs(charge*bmag/mass)
            local kappaPara = n*T*tau/mass
            local kappaPerp = n*T*tau/mass/(Omega*tau)^2
            return kappaPara, kappaPerp
         end
         self._anisotropicDiffusion:setKappaFunction(kappaFunction)
      elseif self._kappaMode=="from-scratch" then
         assert(false,
                self._pfx.."kappaMode 'from-scratch' is not implemented.")
         -- TODO Compute tau and then kappa.
         self._anisotropicDiffusion:setKappaFunction(kappaFunction)
      end

      self._anisotropicDiffusion:setDtAndCflRate(dt)
      self._anisotropicDiffusion:setAuxField(fld)
      self._anisotropicDiffusion:setCalcDivQSwitch(false)
      local myStatus, myDtSuggested = self._anisotropicDiffusion:advance(
         tCurr, {buf, emf, staticEmf}, {buf, buf})
      status = status and myStatus
      dtSuggested = math.min(dtSuggested, myDtSuggested)
      if not status then return status, dtSuggested end

      -- Step 3: Add additional terms to q.
      -- TODO Implement these terms.
      if self._addAdditionalQ then
      end

      -- Step 4: Compute div(q).
      self._anisotropicDiffusion:setCalcQSwitch(false)
      self._anisotropicDiffusion:setCalcDivQSwitch(true)
      self._anisotropicDiffusion:advance(
         tCurr, {buf, emf, staticEmf}, {buf, buf})

      -- Step 5: Add -div(q) onto total energy.
      for idx in localRange:rowMajorIter() do
         fld:fill(fldIdxr(idx), fldPtr)
         buf:fill(bufIdxr(idx), bufPtr)
         fldPtr[5] = fldPtr[5] - bufPtr[5] * dt
      end

   end  -- Loop over species ends.

   return status, dtSuggested
end

function BraginskiiHeatConduction:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiHeatConduction
