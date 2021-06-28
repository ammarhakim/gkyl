-- Gkyl ------------------------------------------------------------------------
--
-- An SSP-RK2 stepper. Mildly unstable and in general should not be used.
-- See gkyl docs for formulas for various SSP-RK schemes:
--   http://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Time = require "Lib.Time"

local SSPRK2 = Proto()

-- Store table passed to it and defer construction.
function SSPRK2:init(tbl)
   self.tbl = tbl

   self.cflFrac = 1.0
   self.numFields = 3
end

function SSPRK2:createSolver(appStatus, stepperFuncs, appsIn)

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
         if stat.dt_actual < dt then
            -- Diagnostics.
            local dt_relDiff  = (dt-stat.dt_actual)/stat.dt_actual
            appStatus.dtDiff[2][1] = math.min(appStatus.dtDiff[2][1], dt_relDiff)
            appStatus.dtDiff[2][2] = math.max(appStatus.dtDiff[2][2], dt_relDiff)
            appStatus.nFail[2]     = appStatus.nFail[2] + 1

            dtNext, nextState = stat.dt_actual, self.RK_STAGE_1
         else
            local tm = Time.clock()
            self.combine(2, 1.0/2.0, 1, 1.0/2.0, 3)
            self.stepperTime = self.stepperTime + (Time.clock() - tm)
            self.copy(1, 2)
            dtNext, nextState = dt, self.RK_COMPLETE
         end
         return dtNext, nextState
      end,
   }

end

function SSPRK2:advance(tCurr, dtIn)

   local dt, rkState = dtIn, self.RK_STAGE_1

   while rkState ~= self.RK_COMPLETE do
      dt, rkState = self.stages[rkState](tCurr, dt)
   end

   return self.stepStatus
end

return SSPRK2
