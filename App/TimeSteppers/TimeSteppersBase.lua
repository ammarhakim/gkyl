-- Gkyl ---------------------------------------
--
-- Base object for the time steppers.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

local TimeSteppersBase = Proto()

function TimeSteppersBase:init(tbl) end
function TimeSteppersBase:createSolver(appStatus, stepperFuncs, appsIn) end
function TimeSteppersBase:advance(tCurr, dtIn) end

-- Assumes the tbl is an array, i.e., all the keys are
-- successive integers - otherwise #tbl will fail.
-- from: https://unendli.ch/posts/2016-07-22-enumerations-in-lua.html
function TimeSteppersBase:enum(tbl)
    local length = #tbl
    local newTbl = {}
    for i = 1, length do
        local v = tbl[i]
        newTbl[v] = i
    end
    return newTbl
end

return TimeSteppersBase
