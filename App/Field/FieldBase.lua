-- Gkyl ------------------------------------------------------------------------
--
-- App support code: Base object for fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell field base classes.
local FieldBase = Proto()
function FieldBase:init(tbl) self.isElliptic = false end
function FieldBase:readRestart() end
function FieldBase:hasEB() return nil, nil end
function FieldBase:setCfl() end
function FieldBase:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function FieldBase:setBasis(basis) self.basis = basis end
function FieldBase:printDevDiagnostics() end
function FieldBase:accumulateCurrent(dt, current, em) end
function FieldBase:applyBcIdx(tCurr, idx) end
function FieldBase:useBoundaryFlux(tCurr, outIdx) end
function FieldBase:suggestDt() return GKYL_MAX_DOUBLE end
function FieldBase:clearCFL() end
function FieldBase:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end
function FieldBase:clearTimers() end


local ExternalFieldBase = Proto()
function ExternalFieldBase:init(tbl) self.isElliptic = false end
function ExternalFieldBase:hasEB() return nil, nil end
function ExternalFieldBase:setCfl() end
function ExternalFieldBase:setIoMethod(ioMethod) self.ioMethod = ioMethod end
function ExternalFieldBase:setBasis(basis) self.basis = basis end
function ExternalFieldBase:readRestart() end
function ExternalFieldBase:printDevDiagnostics() end
function ExternalFieldBase:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end
function ExternalFieldBase:clearTimers() end

-- NoField ---------------------------------------------------------------------
--
-- Represents no field (nothing is evolved or stored).
--------------------------------------------------------------------------------

local NoField = Proto(FieldBase)

-- Methods for no field object.
function NoField:init(tbl) NoField.super.init(tbl) end
function NoField:fullInit(tbl) end
function NoField:hasEB() return nil, nil end
function NoField:setCfl() end
function NoField:setIoMethod(ioMethod) end
function NoField:setBasis(basis) end
function NoField:setGrid(grid) end
function NoField:alloc(nField) end
function NoField:createSolver() end
function NoField:createDiagnostics() end
function NoField:initField() end
function NoField:write(tm) end
function NoField:writeRestart(tm) end
function NoField:readRestart() return 0.0 end
function NoField:rkStepperFields() return {} end
function NoField:suggestDt() end
function NoField:clearCFL() end
function NoField:advance(tCurr, momIn, emIn, emOut) end
function NoField:updateInDirection(dir, tCurr, dt, fIn, fOut)
   return true, GKYL_MAX_DOUBLE
end
function NoField:applyBcIdx(tCurr, idx) end
function NoField:applyBc(tCurr, emIn) end
function NoField:totalSolverTime() return 0.0 end
function NoField:totalBcTime() return 0.0 end
function NoField:energyCalcTime() return 0.0 end
function NoField:copyRk() end
function NoField:combineRk() end
function NoField:printDevDiagnostics() end
function NoField:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end
function NoField:clearTimers() end

return {
   FieldBase         = FieldBase,
   ExternalFieldBase = ExternalFieldBase,
   FuncFieldBase     = ExternalFieldBase,   -- For backwards compatibility.
   NoField           = NoField
}

