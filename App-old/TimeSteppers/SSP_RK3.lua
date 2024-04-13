-- Gkyl ------------------------------------------------------------------------
--
-- An SSP-RK3 stepper.
-- See gkyl docs for formulas for various SSP-RK schemes:
--   http://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local TimeSteppersBase = require "App.TimeSteppers.TimeSteppersBase"

local SSPRK3 = Proto(TimeSteppersBase)

-- Store table passed to it and defer construction.
function SSPRK3:init(tbl)
   self.tbl = tbl

   self.cflFrac = 1.0
   self.numFields = 3
end

function SSPRK3:createSolver(appStatus, stepperFuncs, appsIn)

   self.combine      = stepperFuncs[1]
   self.copy         = stepperFuncs[2]
   self.dydt         = stepperFuncs[3]
   self.forwardEuler = stepperFuncs[4]

   self.stepperTime = 0

   -- Stepper status.
   self.stepStatus = {success   = true,  dt_actual    = 0.,
                      nextState = 0,     dt_suggested = 0.,}

   -- States.
   self.RK_STAGE_1  = 1
   self.RK_STAGE_2  = 2
   self.RK_STAGE_3  = 3
   self.RK_COMPLETE = -1

   local stat = self.stepStatus
   self.stages = {
      [self.RK_STAGE_1] = function(tCurr, dt)
         self.dydt(tCurr, 1, 2)
         self.forwardEuler(tCurr, dt, 1, 2, stat)
         local dtNext, nextState = stat.dt_actual, self.RK_STAGE_2
         return dtNext, nextState
      end,
      [self.RK_STAGE_2] = function(tCurr, dt)
         self.dydt(tCurr+dt, 2, 3)
         self.forwardEuler(tCurr, dt, 2, 3, stat)
         local dtNext, nextState
         -- MF 2021/08/04: Disabling later stage failues for now, so steps are
         --                as it's been in g2 and not quite like it is in g0 now.
         --if stat.dt_actual < dt then
         --   -- Diagnostics.
         --   local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
         --   appStatus.dtDiff[2][1] = math.min(appStatus.dtDiff[2][1], dt_relDiff)
         --   appStatus.dtDiff[2][2] = math.max(appStatus.dtDiff[2][2], dt_relDiff)
         --   appStatus.nFail[2]     = appStatus.nFail[2] + 1

         --   dtNext, nextState = stat.dt_actual, self.RK_STAGE_1
         --else
            local tm = Time.clock()
            self.combine(2, 3.0/4.0, 1, 1.0/4.0, 3)
            self.stepperTime = self.stepperTime + (Time.clock() - tm)
            dtNext, nextState = dt, self.RK_STAGE_3
         --end
         return dtNext, nextState
      end,
      [self.RK_STAGE_3] = function(tCurr, dt)
         self.dydt(tCurr+dt/2, 2, 3)
         self.forwardEuler(tCurr, dt, 2, 3, stat)
         local dtNext, nextState
         -- MF 2021/08/04: Disabling later stage failues for now, so steps are
         --                as it's been in g2 and not quite like it is in g0 now.
         --if stat.dt_actual < dt then
         --   -- Diagnostics.
         --   local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
         --   appStatus.dtDiff[3][1] = math.min(appStatus.dtDiff[3][1], dt_relDiff)
         --   appStatus.dtDiff[3][2] = math.max(appStatus.dtDiff[3][2], dt_relDiff)
         --   appStatus.nFail[3]     = appStatus.nFail[3] + 1

         --   dtNext, nextState = stat.dt_actual, self.RK_STAGE_1
         --else
            local tm = Time.clock()
            self.combine(2, 1.0/3.0, 1, 2.0/3.0, 3)
            self.stepperTime = self.stepperTime + (Time.clock() - tm)
            self.copy(1, 2)
            dtNext, nextState = dt, self.RK_COMPLETE
         --end
         return dtNext, nextState
      end,
   }

end

function SSPRK3:advance(tCurr, dtIn)

   local dt, rkState = dtIn, self.RK_STAGE_1

   while rkState ~= self.RK_COMPLETE do
      dt, rkState = self.stages[rkState](tCurr, dt)
   end

   return self.stepStatus
end

return SSPRK3
