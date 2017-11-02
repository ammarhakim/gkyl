-- Gkyl ------------------------------------------------------------------------
--
-- Helper functions for use in updaters.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Time = require "Lib.Time"
local Mpi = require "Comm.Mpi"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local _M = {}

function _M.setup(self, tbl)
   assert(tbl.onGrid, "Updater.Base: Must provide grid object using 'onGrid'")
   self._comm = tbl.onGrid:commSet().comm
   self._sharedComm = tbl.onGrid:commSet().sharedComm
   self._nodeComm = tbl.onGrid:commSet().nodeComm
end

-- This function wraps an updater's advance() function and adds a
-- "totalTime" field to the updater, and also synchronizes the status
-- and time-step suggestion across processors.
function _M.advanceFuncWrap(advanceFunc)
   return function(self, tCurr, dt, inFld, outFld)
      self.totalTime = self.totalTime and self.totalTime or 0.0
      local tmStart = Time.clock()
      -- >>> Take the time-step
      local _status, _dtSuggested = advanceFunc(self, tCurr, dt, inFld, outFld)
      -- >>> Done with advance
      self.totalTime = self.totalTime + (Time.clock()-tmStart)

      -- reduce across processors ...
      local myStatus, myDtSuggested = new("int[1]"), new("double[1]")
      local status, dtSuggested = new("int[1]"), new("double[1]")      
      myStatus[0] = _status and 1 or 0; myDtSuggested[0] = _dtSuggested

      Mpi.Allreduce(myStatus, status, 1, Mpi.INT, Mpi.MIN, self._nodeComm)
      Mpi.Allreduce(myDtSuggested, dtSuggested, 1, Mpi.DOUBLE, Mpi.MIN, self._nodeComm)

      return status[0] == 1 and true or false, dtSuggested[0]
   end
end

return _M


