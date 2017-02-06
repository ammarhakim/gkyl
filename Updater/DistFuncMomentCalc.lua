-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Moments updater object
local DistFuncMomentCalc = {}

function DistFuncMomentCalc:new(tbl)
   local self = setmetatable({}, DistFuncMomentCalc)
   Base.setup(self, tbl) -- setup base object

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(DistFuncMomentCalc, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   
   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
DistFuncMomentCalc.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   DistFuncMomentCalc = DistFuncMomentCalc
}
