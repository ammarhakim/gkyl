-- Gkyl ------------------------------------------------------------------------
--
-- Helper functions for use in updaters.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Time = require "Lib.Time"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

_M = {}

-- This function wraps an updater's advance() function and addes a
-- "totalTime" field to the updater, and also synchronizes the status
-- and time-step suggestion across processors.
function _M.advanceFuncWrap(advanceFunc)
   return function(self, tCurr, dt, inFld, outFld)
      self.totalTime = self.totalTime and self.totalTime or 0.0

      local tmStart = Time.clock()
      local status, dtSuggested = advanceFunc(self, tCurr, dt, inFld, outFld)
      self.totalTime = self.totalTime + (Time.clock()-tmStart)

      return status, dtSuggested
   end
end

return _M


