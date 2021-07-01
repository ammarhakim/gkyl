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
   self.numFields = 4
   self.numStates = 5
end

-- This set of functions determines factors which feed into RK scheme
-- (Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal of
-- Computational Physics, 257(PA), 594–626. doi:10.1016/j.jcp.2013.08.021).
local function mu(j)
   return (2*j-1)/j   -- RKL1.
--   return (2*j-1)*(j+2)*((j-1)^2)/(j*(j-2)*((j+1)^2))   -- RKL2.
end
local function nu(j)
   return (1-j)/j   -- RKL1.
--   return -((j-1)^3)*(j^2-4)/((j^3)*(j-1)*(j+3))   -- RKL2.
end
local function muTilde(s,j)
   -- RKL1.
   return ((2*j-1)/j)*(2/(s^2+s))
   -- RKL2.
--   return (4/(3*(s^2+s-2)))*(j==1 and 1/3 or ((2*j-1)*(j+2)*((j-1)^2)/(j*(j-2)*((j+1)^2)))*(2/(s^2+s)))
end
local function gammaTilde(j)   -- For RKL2.
   return -(j-1)*(j+2)*(2*j-1)*(j^2-j+2)/(2*(j^2)*(j-2)*((j+1)^2))
end

local function calcNumStages(dtRatio)
   -- dtRatio is the ratio of the time step by which we wish to step
   -- the split (e.g. parabolic) operator to the time step of the non-split
   -- (e.g. hyperbolic) operators.
   -- Following Meyer et al. MNRAS 422, 2102–2115 (2012) we choose the smallest odd
   -- integer that satisfies the inequality in equation 26 of that paper.
   -- RKL1.
   --return math.ceil((1/2)*(math.sqrt(1+8*dtRatio)-1))
   return math.floor(math.ceil(1/2*(math.sqrt(1+8*dtRatio)-1))/2)*2+1
   -- RKL2.
   --return math.ceil((1/2)*(math.sqrt(9+16*dtRatio)-1))
end

local function dtMaxAllowed(s, dtExplicit)
   -- RKL1, eq. 11 of Meyer et al, JCP (2014).
   return dtExplicit*(s^2+s)/2
   -- RKL2, eq. 19 of Meyer et al, JCP (2014).
   --return dtExplicit*(s^2+s-2)/4
end

local function bufferIdxs(idxIn)
   -- Given a field index (e.g. in species' :rkStepperFields), return
   -- the indices of the other fields which we may use as temporary buffers.
   if idxIn==1 then     return {2,3,4}
   elseif idxIn==2 then return {3,4,1}
   elseif idxIn==3 then return {4,1,2}
   elseif idxIn==4 then return {1,2,3} end
end

function OperatorSplitSSPRK3:sts(tCurr, outIdx, dtIn, inIdx, stat)
   -- Takes f[inIdx] and computes f[outIdx]=f[inIdx]+dt*dydtSplit with super time
   -- stepping, where dydtSplit is the time rate of change due to the split operators.
   -- IMPORTANT: Below there are some /2 or *2. because we assume sts is called with dt/2.

   local fIdxs = bufferIdxs(inIdx)
   local jm1, jm2, fDotIdx = fIdxs[1], fIdxs[2], fIdxs[3]

   -- Stage 1. See equation 9 of Meyer et al. JCP 2014.
   self.dydtSplit(tCurr, inIdx, fDotIdx)

   local dt, dtSplitMin = dtIn, GKYL_MAX_DOUBLE
   for _, s in pairs(self.species) do dtSplitMin = math.min(dtSplitMin, s:suggestDt()) end

   local numStages = calcNumStages(dt/dtSplitMin)   -- Number of RKL stages.

   local dtMax = dtMaxAllowed(numStages, dtSplitMin)
   if dtMax < 0.99*dt then dt = dtMax end

   stat.dt_actual    = dt*2.
   stat.dt_suggested = dt*2.

   for _, s in lume.orderedIter(self.species) do
      local fIn = s:rkStepperFields()[inIdx]
      local fDot = s:rkStepperFields()[fDotIdx]
      local fjm1, fjm2 = s:rkStepperFields()[jm1], s:rkStepperFields()[jm2]
      fjm2:copy(fIn)
      -- The following can rewrite f[jm1] so muTilde*dt*fDot has to go first in :combine.
      fjm1:combine(muTilde(numStages,1)*dt, fDot, 1.0, fIn)
   end
   for _, s in lume.orderedIter(self.species) do
      s:applyBcIdx(tCurr, self.field, self.externalField, inIdx, jm1)
   end

   for j = 2, numStages do   -- Remaining stages.
      self.dydtSplit(tCurr, jm1, fDotIdx)

      for _, s in lume.orderedIter(self.species) do
         local fDot = s:rkStepperFields()[fDotIdx]
         local fjm1, fjm2 = s:rkStepperFields()[jm1], s:rkStepperFields()[jm2]
         local fj = fDot
         
         -- The following rewrites f[fDotIdx], so muTilde*dt*fDot has to go first in :combine.
         fj:combine(muTilde(numStages,j)*dt, fDot, mu(j), fjm1, nu(j), fjm2)
      end
      for _, s in lume.orderedIter(self.species) do
         s:applyBcIdx(tCurr, self.field, self.externalField, jm1, fDotIdx)
      end

      -- Reset fields for next stage.
      for _, s in lume.orderedIter(self.species) do
         local fj, fjm1, fjm2 = s:rkStepperFields()[fDotIdx], s:rkStepperFields()[jm1], s:rkStepperFields()[jm2]
         fjm2:copy(fjm1)
         fjm1:copy(fj)
      end
   end

   for _, s in lume.orderedIter(self.species) do
      local fj, fOut = s:rkStepperFields()[jm1], s:rkStepperFields()[outIdx]
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
--         local dtNext, nextState = dt, self.RK_STAGE_1
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
--   print(string.format("pre dt=%g | dtSplitMin=%g | numStages=%d | dt_act=%g",dt, dtSplitMin, numStages, stat.dt_actual))
--   print(string.format("aft dt=%g | dtSplitMin=%g | numStages=%d | dtMax=%g | dt_act=%g",dt, dtSplitMin, numStages, dtMax, stat.dt_actual))
--   print(" ")
