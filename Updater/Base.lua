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
   self._myStatus, self._myDtSuggested = ffi.new("int[2]"), ffi.new("double[2]")
   self._status, self._dtSuggested = ffi.new("int[2]"), ffi.new("double[2]")
end

-- must be provided by derived objects
function _M:_advance(tCurr, dt, inFld, outFld)
   assert(true, "_advance method not provided!")
end

-- return various comms for communications
function _M:getNodeComm()
   return self._nodeComm
end
function _M:getComm()
   return self._comm
end
function _M:getSharedComm()
   return self._sharedComm
end

-- This function wraps derived updater's _advance() function and
-- computes a "totalTime", and also synchronizes the status and
-- time-step suggestion across processors.
function _M:advance(tCurr, dt, inFld, outFld)

   -- Take the time-step, measuring how long it took
   local tmStart = Time.clock()
   local status, dtSuggested = self:_advance(tCurr, dt, inFld, outFld)
   self.totalTime = self.totalTime + (Time.clock()-tmStart)

   -- reduce across processors ...
   self._myStatus[0] = status and 1 or 0
   self._myDtSuggested[0] = dtSuggested

   Mpi.Allreduce(self._myStatus, self._status, 1, Mpi.INT, Mpi.LAND, self._comm)
   Mpi.Allreduce(self._myDtSuggested, self._dtSuggested, 1, Mpi.DOUBLE, Mpi.MIN, self._comm)
   
   return self._status[0] == 1 and true or false, self._dtSuggested[0]
end

return _M


