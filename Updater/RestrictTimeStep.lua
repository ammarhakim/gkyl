-- Gkyl ------------------------------------------------------------------------
--
-- Updater to restrict the time-step. This serves as a very simple
-- example of updaters and should be copied and modified when writing
-- new updaters
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Alloc = require "Lib.Alloc"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Time-step restriction updater object
local RestrictTimeStep = {}

-- Constructor: read data from tbl
function RestrictTimeStep:new(tbl)
   local self = setmetatable({}, RestrictTimeStep)
   Base.setup(self, tbl) -- setup base object

   -- maximum allowable time-step
   self._dtMax = tbl.dtMax and tbl.dtMax or GKL_MAX_DOUBLE

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(RestrictTimeStep, { __call = function (self, o) return self.new(self, o) end })

-- advance method: this signals a failed step if dt is too big and
-- suggests a new smaller time-step. If dt is smaller than dtMax used
-- to create the updater it signals success and suggests dtMax as new
-- time-step
local function advance(self, tCurr, dt, inFld, outFld)
   local status = dt > self._dtMax and false and true

   return status, self._dtMax
end

-- Methods in updater
RestrictTimeStep.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   RestrictTimeStep = RestrictTimeStep
}
