-- Gkyl ------------------------------------------------------------------------
--
-- Helper functions for use in updaters.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Mpi   = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time  = require "Lib.Time"
local ffi   = require "ffi"
local xsys  = require "xsys"

local _M = Proto()

function _M:init(tbl)
   -- Some updaters use the communicators to perform parallel operations.
   local onGrid = tbl.onGrid
   if onGrid then
      self._comm = onGrid:commSet().comm
   end

   self.totalTime = 0.0  -- Time taken by the advance method.

   -- Some (FV) updaters are used for time stepping and require parallel
   -- reduction of stepping status and size.
   self._myStatus, self._myDtSuggested = ffi.new("int[2]"), ffi.new("double[2]")
   self._status, self._dtSuggested     = ffi.new("int[2]"), ffi.new("double[2]")

   self._dt = 0.0
   self._cflRateByCell = nil

   -- check if we should skip the device update for this specific
   -- updater
   local skipDevice = xsys.pickBool(tbl.skipDevice, false)

   -- store advance method to use
   self._advanceFunc = self._advance
   if GKYL_USE_GPU then
      if self._advanceOnDevice then
	 self._advanceFunc = self._advanceOnDevice
      else
	 -- this means we don't really have a GPU implementation of
	 -- the updater: needs to be handled carefully
	 self._advanceFunc = self._advanceNoDeviceImpl
      end
   end
end

function _M:_advanceNoDeviceImpl(tCurr, inFld, outFld)
   -- Copy input fields from device -> host.
   for _, fld in ipairs(inFld) do 
      if type(fld)=="table" and fld._zero then fld:copyDeviceToHost() end
   end
    -- Also copy output fields in case they are inputs too,
    -- or are incremented rather than overwritten.
   for _, fld in ipairs(outFld) do
      if type(fld)=="table" and fld._zero then fld:copyDeviceToHost() end
   end

   self:_advance(tCurr, inFld, outFld)

   -- Copy output fields from host -> device.
   for _, fld in ipairs(outFld) do 
      if type(fld)=="table" and fld._zero then fld:copyHostToDevice() end
   end
   -- Also copy input fields in case they were modified.
   for _, fld in ipairs(inFld) do 
      if type(fld)=="table" and fld._zero then fld:copyHostToDevice() end
   end
end

-- must be provided by derived objects
function _M:_advance(tCurr, inFld, outFld)
   assert(true, "_advance method not provided!")
end

-- return various comms for communications
function _M:getComm() return self._comm end

-- This function wraps derived updater's _advance() function and
-- computes a "total Time", and also synchronizes the status and
-- time-step suggestion across processors.
function _M:advance(tCurr, inFld, outFld)
   local tmStart = Time.clock()

   local status, dtSuggested = self:_advanceFunc(tCurr, inFld, outFld)

   self.totalTime = self.totalTime + (Time.clock()-tmStart)
   return status, dtSuggested
end

function _M:reduceStatusDt(status, dtSuggested)
   if status ~= nil and dtSuggested ~= nil then 
      self._myStatus[0]      = status and 1 or 0
      self._myDtSuggested[0] = dtSuggested

      Mpi.Allreduce(self._myStatus, self._status, 1, Mpi.INT, Mpi.LAND, self._comm)
      Mpi.Allreduce(self._myDtSuggested, self._dtSuggested, 1, Mpi.DOUBLE, Mpi.MIN, self._comm)
      
      return self._status[0] == 1 and true or false, self._dtSuggested[0]
   end
end

-- set up pointers to dt and cflRateByCell
function _M:setDtAndCflRate(dt, cflRateByCell)
   self._dt = dt
   self._cflRateByCell = cflRateByCell
end

function _M:printDevDiagnostics() end

return _M


