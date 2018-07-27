-- Gkyl ------------------------------------------------------------------------
--
-- App support code: Base object for fields
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- empty shell field base classes
local FieldBase = Proto()
function FieldBase:init(tbl)
   self.isElliptic = false
end

local FuncFieldBase = Proto()
function FuncFieldBase:init(tbl)
   self.isElliptic = false
end

-- NoField ---------------------------------------------------------------------
--
-- Represents no field (nothing is evolved or stored)
--------------------------------------------------------------------------------

local NoField = Proto(FieldBase)

-- methods for no field object
function NoField:init(tbl)
   NoField.super.init(tbl)
end
function NoField:fullInit(tbl) end
function NoField:hasEB() return false, false end
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
function NoField:forwardEuler(tCurr, dt, momIn, emIn, emOut)
   return true, GKYL_MAX_DOUBLE
end
function NoField:updateInDirection(dir, tCurr, dt, fIn, fOut)
   return true, GKYL_MAX_DOUBLE
end
function NoField:applyBc(tCurr, dt, emIn) end
function NoField:totalSolverTime() return 0.0 end
function NoField:totalBcTime() return 0.0 end
function NoField:energyCalcTime() return 0.0 end
function NoField:copyRk() end
function NoField:combineRk() end

return {
   FieldBase = FieldBase,
   FuncFieldBase = FuncFieldBase,
   NoField = NoField
}

