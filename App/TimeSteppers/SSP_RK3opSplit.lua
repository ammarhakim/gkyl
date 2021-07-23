-- Gkyl ------------------------------------------------------------------------
--
-- An SSP-RK3 stepper with an operator (Strang) splitting, applying
-- a dt/2 update of certain terms before an after the SSP-RK3.
-- See gkyl docs for formulas for various SSP-RK schemes:
--   http://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local lume = require "Lib.lume"

local OperatorSplitSSPRK3 = Proto()

-- Store table passed to it and defer construction.
function OperatorSplitSSPRK3:init(tbl)
   self.tbl = tbl

   self.cflFrac = 1.0
   self.numStates = 5

   self.isRKL1 = false   -- True: RKL1. False: RKL2.

   self.numFields = self.isRKL1 and 4 or 5
end

-- This set of functions determines factors which feed into RK scheme
-- (Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal of
-- Computational Physics, 257(PA), 594–626. doi:10.1016/j.jcp.2013.08.021).
local function bRKL2(j, isRKL1)   -- For RKL2.
   return j<3 and 1/3 or (j^2+j-2)/(2*j*(j+1))
end
local function mu(j, isRKL1)
   return isRKL1
      and (2*j-1)/j   -- RKL1.
      or ((2*j-1)/j)*(bRKL2(j)/bRKL2(j-1))   -- RKL2.
end
local function nu(j,isRKL1)
   return isRKL1
      and (1-j)/j   -- RKL1.
      or -((j-1)/j)*(bRKL2(j)/bRKL2(j-2))   -- RKL2.
end
local function muTilde(s,j,isRKL1)
   return isRKL1
      and ((2*j-1)/j)*(2/(s^2+s))   -- RKL1.
      or (4/(s^2+s-2))*(j==1 and 1/3 or ((2*j-1)/j)*(bRKL2(j)/bRKL2(j-1)))   -- RKL2.
end
local function gammaTilde(s,j)   -- For RKL2.
   return -(1-bRKL2(j-1))*muTilde(s,j,false)
end

local function calcNumStages(dtRatio, isRKL1)
   -- dtRatio is the ratio of the time step by which we wish to step
   -- the split (e.g. parabolic) operator to the time step of the non-split
   -- (e.g. hyperbolic) operators.
   -- Following Meyer et al. MNRAS 422, 2102–2115 (2012) we choose the smallest odd
   -- integer that satisfies the inequality in equation 26 of that paper.
   return isRKL1
      and math.floor(math.ceil(0.5*(math.sqrt(1+8*dtRatio)-1))/2)*2+1   -- RKL1.
      or math.floor(math.ceil(0.5*(math.sqrt(9+16*dtRatio)-1))/2)*2+1   -- RKL2.
end

local function dtMaxAllowed(s, dtExplicit, isRKL1)
   return isRKL1
     and dtExplicit*(s^2+s)/2   -- RKL1, eq. 11 of Meyer et al, JCP (2014).
      or dtExplicit*(s^2+s-2)/4   -- RKL2, eq. 19 of Meyer et al, JCP (2014).
end

local function bufferIdxs(idxIn, isRKL1)
   -- Given a field index (e.g. in species' :rkStepperFields), return
   -- the indices of the other fields which we may use as temporary buffers.
   if isRKL1 then
      if idxIn==1 then     return {2,3,4}
      elseif idxIn==2 then return {3,4,1}
      elseif idxIn==3 then return {4,1,2}
      elseif idxIn==4 then return {1,2,3} end
   else
      if idxIn==1 then     return {2,3,4,5}
      elseif idxIn==2 then return {3,4,5,1}
      elseif idxIn==3 then return {4,5,1,2}
      elseif idxIn==4 then return {5,1,2,3}
      elseif idxIn==5 then return {1,2,3,4} end
   end
end

function OperatorSplitSSPRK3:sts(tCurr, outIdx, dtIn, inIdx, stat)
   -- Takes f[inIdx] and computes f[outIdx]=f[inIdx]+dt*dydtSplit with super time
   -- stepping, where dydtSplit is the time rate of change due to the split operators.

   local isRKL1 = self.isRKL1

   local fIdxs = bufferIdxs(inIdx, isRKL1)
   local jm1, jm2, fDotIdx, fDot0Idx = fIdxs[1], fIdxs[2], fIdxs[3], fIdxs[4]

   -- Stage 1. See equation 9 of Meyer et al. JCP 2014.
   self.dydtSplit(tCurr, inIdx, fDotIdx)

   -- Below :suggestDtSplit sets the largest dt by which sts should step the solution.
   local dt, dtSplitExp = dtIn, dtIn
   for _, s in pairs(self.species) do
      dt = math.min(dt, s:suggestDtSplit())
      dtSplitExp = math.min(dtSplitExp, s:suggestDt())
   end
   -- IMPORTANT: *2 below because we assume sts is called with dt/2.
   stat.dt_actual, stat.dt_suggested = dt*2, dt*2

   local numStages = calcNumStages(dt/dtSplitExp, isRKL1)   -- Number of RKL stages.
--   print(string.format("dt=%g | dtSplitExp=%g | numStages=%d",dt,dtSplitExp,numStages))

   for _, s in lume.orderedIter(self.species) do
      local flds = s:rkStepperFields()
      local fIn, fDot, fjm1, fjm2 = flds[inIdx], flds[fDotIdx], flds[jm1], flds[jm2]
      fjm2:copy(fIn)
      fjm1:combine(muTilde(numStages,1,isRKL1)*dt, fDot, 1.0, fIn)

      if not isRKL1 then   -- For RKL2.
         local fDot0 = s:rkStepperFields()[fDot0Idx]
         fDot0:copy(fDot)
      end
   end
   for _, s in lume.orderedIter(self.species) do
      s:applyBcIdx(tCurr, self.field, self.externalField, inIdx, jm1)
   end

   for jS = 2, numStages do   -- Remaining stages.
      self.dydtSplit(tCurr, jm1, fDotIdx)

      for _, s in lume.orderedIter(self.species) do
         local flds = s:rkStepperFields()
         local fDot, fjm1, fjm2 = flds[fDotIdx], flds[jm1], flds[jm2]
         local fj = fDot
         -- The following rewrites f[fDotIdx], so muTilde*dt*fDot has to go first in :combine.
         fj:combine(muTilde(numStages,jS,isRKL1)*dt, fDot, mu(jS,isRKL1), fjm1, nu(jS,isRKL1), fjm2)
            
         if not isRKL1 then   -- For RKL2.
            local fDot0 = s:rkStepperFields()[fDot0Idx]
            local fIn = s:rkStepperFields()[inIdx]
            fj:accumulate(1.-mu(jS,isRKL1)-nu(jS,isRKL1), fIn, gammaTilde(numStages,jS)*dt, fDot0)
         end
      end
      for _, s in lume.orderedIter(self.species) do
         s:applyBcIdx(tCurr, self.field, self.externalField, jm1, fDotIdx)
      end

      -- Reset fields for next stage.
      for _, s in lume.orderedIter(self.species) do
         local flds = s:rkStepperFields()
         local fj, fjm1, fjm2 = flds[fDotIdx], flds[jm1], flds[jm2]
         fjm2:copy(fjm1)
         fjm1:copy(fj)
      end
   end

   for _, s in lume.orderedIter(self.species) do
      local flds = s:rkStepperFields()
      local fj, fOut = flds[jm1], flds[outIdx]
      fOut:copy(fj)
   end
end

function OperatorSplitSSPRK3:createSolver(appStatus, stepperFuncs, appsIn)

   self.combine      = stepperFuncs[1]
   self.copy         = stepperFuncs[2]
   self.dydt         = stepperFuncs[3]
   self.forwardEuler = stepperFuncs[4]
   self.dydtSplit    = stepperFuncs[5]

   self.species       = appsIn[1]
   self.field         = appsIn[2]
   self.externalField = appsIn[3]

   self.stepperTime = 0

   -- Stepper status.
   self.stepStatus = {success   = true,  dt_actual    = 0.,
                      nextState = 0,     dt_suggested = 0.,}

   -- States.
   self.SPLIT_STAGE_1  = 0
   self.RK_STAGE_1     = 1
   self.RK_STAGE_2     = 2
   self.RK_STAGE_3     = 3
   self.SPLIT_STAGE_2  = 4
   self.RK_COMPLETE    = -1

   local stat = self.stepStatus
   self.stages = {
      [self.SPLIT_STAGE_1] = function(tCurr, dt)
         self:sts(tCurr, 1, dt/2., 1, stat)
         local dtNext, nextState = stat.dt_actual, self.RK_STAGE_1
         return dtNext, nextState
      end,
      [self.RK_STAGE_1] = function(tCurr, dt)
         self.dydt(tCurr, 1, 2)
         self.forwardEuler(tCurr, dt, 1, 2, stat)
         local dtNext, nextState
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff = (dt-stat.dt_actual)/stat.dt_actual
            appStatus.dtDiff[2][1] = math.min(appStatus.dtDiff[2][1], dt_relDiff)
            appStatus.dtDiff[2][2] = math.max(appStatus.dtDiff[2][2], dt_relDiff)
            appStatus.nFail[2]     = appStatus.nFail[2] + 1

            dtNext, nextState = stat.dt_actual, self.SPLIT_STAGE_1
         else
            dtNext, nextState = stat.dt_actual, self.RK_STAGE_2
         end
         return dtNext, nextState
      end,
      [self.RK_STAGE_2] = function(tCurr, dt)
         self.dydt(tCurr+dt, 2, 3)
         self.forwardEuler(tCurr, dt, 2, 3, stat)
         local dtNext, nextState
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff = (dt-stat.dt_actual)/stat.dt_actual
            appStatus.dtDiff[3][1] = math.min(appStatus.dtDiff[3][1], dt_relDiff)
            appStatus.dtDiff[3][2] = math.max(appStatus.dtDiff[3][2], dt_relDiff)
            appStatus.nFail[3]     = appStatus.nFail[3] + 1

            dtNext, nextState = stat.dt_actual, self.SPLIT_STAGE_1
         else
            local tm = Time.clock()
            self.combine(2, 3.0/4.0, 1, 1.0/4.0, 3)
            self.stepperTime = self.stepperTime + (Time.clock() - tm)
            dtNext, nextState = dt, self.RK_STAGE_3
         end
         return dtNext, nextState
      end,
      [self.RK_STAGE_3] = function(tCurr, dt)
         self.dydt(tCurr+dt/2, 2, 3)
         self.forwardEuler(tCurr, dt, 2, 3, stat)
         local dtNext, nextState
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff = (dt-stat.dt_actual)/stat.dt_actual
            appStatus.dtDiff[4][1] = math.min(appStatus.dtDiff[4][1], dt_relDiff)
            appStatus.dtDiff[4][2] = math.max(appStatus.dtDiff[4][2], dt_relDiff)
            appStatus.nFail[4]     = appStatus.nFail[4] + 1

            dtNext, nextState = stat.dt_actual, self.SPLIT_STAGE_1
         else
            local tm = Time.clock()
            self.combine(2, 1.0/3.0, 1, 2.0/3.0, 3)
            self.stepperTime = self.stepperTime + (Time.clock() - tm)
            self.copy(1, 2)
            dtNext, nextState = dt, self.SPLIT_STAGE_2
         end
         return dtNext, nextState
      end,
      [self.SPLIT_STAGE_2] = function(tCurr, dt)
         self:sts(tCurr, 1, dt/2., 1, stat)
         local dtNext, nextState
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff = (dt-stat.dt_actual)/stat.dt_actual
            appStatus.dtDiff[5][1] = math.min(appStatus.dtDiff[5][1], dt_relDiff)
            appStatus.dtDiff[5][2] = math.max(appStatus.dtDiff[5][2], dt_relDiff)
            appStatus.nFail[5]     = appStatus.nFail[5] + 1

            dtNext, nextState = stat.dt_actual, self.SPLIT_STAGE_1
         else
            dtNext, nextState = dt, self.RK_COMPLETE
         end
         local dtNext, nextState = dt, self.RK_COMPLETE
         return dtNext, nextState
      end,
   }

end

function OperatorSplitSSPRK3:advance(tCurr, dtIn)

   local dt, rkState = dtIn, self.SPLIT_STAGE_1

   while rkState ~= self.RK_COMPLETE do
      dt, rkState = self.stages[rkState](tCurr, dt)
   end

   return self.stepStatus
end

return OperatorSplitSSPRK3
