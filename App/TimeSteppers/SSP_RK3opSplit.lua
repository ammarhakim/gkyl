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
   self.numFields = 3
   self.numStates = 5
end

-- This set of functions determines factors which feed into RK scheme
-- (Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal of
-- Computational Physics, 257(PA), 594–626. doi:10.1016/j.jcp.2013.08.021).
local function mu(j) return (2*j-1)/j end
local function nu(j) return (1-j)/j end
local function mubar(s,j) return (2*j-1)/j*2/(s^2+s) end

local function calcNumStages(dtRatio)
   -- dtRatio is the ratio of the time step by which we wish to step
   -- the split (e.g. parabolic) operator to the time step of the non-split
   -- (e.g. hyperbolic) operators.
   -- RKL1.
   return math.ceil(1/2*(math.sqrt(1+8*dtRatio)-1))
   -- RKL2.
   --return math.ceil(1/2*(math.sqrt(9+16*dtRatio)-1))
end

local function extraFields(idxIn)
   -- Given a field index (e.g. in species' :rkStepperFields), return
   -- the indices of the other fields which we may use as temporary
   -- buffers.
   if idxIn==1 then     return {2,3}
   elseif idxIn==2 then return {3,1}
   elseif idxIn==3 then return {1,2} end
end

function OperatorSplitSSPRK3:sts(tCurr, outIdx, dt, fDotIdx, inIdx)
   -- Takes f[inIdx] and fDot=df/dt=f[fDotIdx] and computes f[outIdx]=f[inIdx]+dt*fDot.
   -- That is, it takes a forward euler step with super time stepping.

   -- Determine the time step supported by the split operators. 
   local dtSplitMin = GKYL_MAX_DOUBLE
   for _, s in pairs(self.species) do dtSplitMin = math.min(dtSplitMin, s:suggestDt()) end

   local numStages = calcNumStages(dt/dtSplitMin)   -- Number of RKL stages.

   local idx0, idx1, idx2 = inIdx, extraFields(inIdx)[1], extraFields(inIdx)[2]

   -- Stage 1 (initial fDot is already computed).
   for _, s in lume.orderedIter(self.species) do
      local fIn = s:rkStepperFields()[inIdx]
      local fDot0 = s:rkStepperFields()[fDotIdx]
      local fJ1, fJ2 = s:rkStepperFields()[idx1], s:rkStepperFields()[idx2]
      fJ2:copy(fIn)
      fJ1:combine(1.0, fIn, mubar(numStages,1)*dt, fDot0)
   end
   for _, s in lume.orderedIter(self.species) do
      s:applyBcIdx(tCurr, self.field, self.externalField, idx1, idx1)
   end

   -- Remaining stages.
   for j = 2, numStages do
      self.dydtSplit(tCurr, idx1, idx0)

      for _, s in lume.orderedIter(self.species) do
         local fDot = s:rkStepperFields()[idx0]
         local fJ, fJ1, fJ2 = fDot, s:rkStepperFields()[idx1], s:rkStepperFields()[idx2]
         fJ:combine(mu(j), fJ1, nu(j), fJ2, mubar(numStages,j)*dt, fDot)
      end

      for _, s in lume.orderedIter(self.species) do
         s:applyBcIdx(tCurr, self.field, self.externalField, idx0, idx0)
      end
      -- Reset fields for next stage.
      for _, s in lume.orderedIter(self.species) do
         local fJ, fJ1, fJ2 = s:rkStepperFields()[idx0], s:rkStepperFields()[idx1], s:rkStepperFields()[idx2]
         fJ2:copy(fJ1)
         fJ1:copy(fJ)
      end
   end

   for _, s in lume.orderedIter(self.species) do
      local fJ = s:rkStepperFields()[idx0]
      local fOut = s:rkStepperFields()[outIdx]
      fOut:copy(fJ)
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
         self.dydtSplit(tCurr, 1, 2)
         self:sts(tCurr, 1, dt/2., 2, 1)
         local dtNext, nextState = dt, self.RK_STAGE_1
         return dtNext, nextState
      end,
      [self.RK_STAGE_1] = function(tCurr, dt)
         self.dydt(tCurr, 1, 2)
         self.forwardEuler(tCurr, dt, 1, 2, stat)
         local dtNext, nextState
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
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
            local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
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
            local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
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
         self.dydtSplit(tCurr, 1, 2)
         self:sts(tCurr, 1, dt/2., 2, 1)
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
