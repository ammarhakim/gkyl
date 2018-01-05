-- Gkyl ------------------------------------------------------------------------
--
-- Helper functions for use in updaters.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"

-- system libraries
local ffi = require "ffi"

local _M = Proto()

function _M:init(tbl)
   assert(tbl.onGrid, "Updater.Base: Must provide grid object using 'onGrid'")
   self._comm = tbl.onGrid:commSet().comm
   self._sharedComm = tbl.onGrid:commSet().sharedComm
   self._nodeComm = tbl.onGrid:commSet().nodeComm
   self.totalTime = 0.0 
end

-- must be provided by derived objects
function _M:_advance(tCurr, dt, inFld, outFld)
   assert(true, "_advance method not provided!")
end

-- This function wraps derived updater's _advance() function and
-- computes a "totalTime", and also synchronizes the status and
-- time-step suggestion across processors.
function _M:advance(tCurr, dt, inFld, outFld)
   local tmStart = Time.clock()
   -- >>> Take the time-step
   local _status, _dtSuggested = self:_advance(tCurr, dt, inFld, outFld)
   -- >>> Done with advance
   self.totalTime = self.totalTime + (Time.clock()-tmStart)

   -- reduce across processors ...
   local myStatus, myDtSuggested = ffi.new("int[1]"), ffi.new("double[1]")
   local status, dtSuggested = ffi.new("int[1]"), ffi.new("double[1]")      
   myStatus[0] = _status and 1 or 0; myDtSuggested[0] = _dtSuggested
   
   Mpi.Allreduce(myStatus, status, 1, Mpi.INT, Mpi.MIN, self._nodeComm)
   Mpi.Allreduce(myDtSuggested, dtSuggested, 1, Mpi.DOUBLE, Mpi.MIN, self._nodeComm)
   
   return status[0] == 1 and true or false, dtSuggested[0]
end

return _M


