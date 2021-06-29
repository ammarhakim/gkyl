-- Gkyl ------------------------------------------------------------------------
--
-- An RK1 stepper. Unstable! Only for testing.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Time = require "Lib.Time"

local RK1 = Proto()

-- Store table passed to it and defer construction.
function RK1:init(tbl)
   self.tbl = tbl

   self.cflFrac = 1.0
   self.numFields = 3
   self.numStates = 1
end

function RK1:createSolver(appStatus, stepperFuncs, appsIn)

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
   self.RK_COMPLETE = -1

   local stat = self.stepStatus
   self.stages = {
      [self.RK_STAGE_1] = function(tCurr, dt)
         self.dydt(tCurr, 1, 2)
         self.forwardEuler(tCurr, dt, 1, 2, stat)
         self.copy(1, 2)
         local dtNext, nextState = stat.dt_actual, self.RK_COMPLETE
         return dtNext, nextState
      end,
   }

end

function RK1:advance(tCurr, dtIn)

   local dt, rkState = dtIn, self.RK_STAGE_1

   while rkState ~= self.RK_COMPLETE do
      dt, rkState = self.stages[rkState](tCurr, dt)
   end

   return self.stepStatus
end

return RK1
