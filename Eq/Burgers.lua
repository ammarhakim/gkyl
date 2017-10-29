-- Gkyl ------------------------------------------------------------------------
--
-- Burgers equation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Burgers equation
local Burgers = {}

-- constructor
function Burgers:new(tbl)
   local self = setmetatable({}, Burgers)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Burgers, { __call = function (self, o) return self.new(self, o) end })

-- Methods 
Burgers.__index = {
   numEquations = function (self) return 1 end,
   numWaves = function (self) return 1 end,
   flux = function (self, dir, qIn, fOut)
      fOut[1] = 0.5*qIn[1]*qIn[1]
   end,
   isPositive = function (self, q)
      return true
   end,
   rp = function (self, dir, delta, ql, qr, waves, s)
      waves[1][1] = delta[1]
      s[1] = 0.5*(ql[1]+qr[1])
   end,
   qFluctuations = function (self, dir, ql, qr, waves, s, amdq, apdq)
      amdq[1], apdq[1] = 0.0, 0.0
      if s[1] < 0 then
	 amdq[1] = s[1]*waves[1][1]
      else
	 apdq[1] = s[1]*waves[1][1]
      end
   end,
}

return {
   Burgers = Burgers,
}


