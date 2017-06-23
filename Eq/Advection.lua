-- Gkyl ------------------------------------------------------------------------
--
-- Linear, constant velocity advection equation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Constant velocity advection equation in 
local Advection = {}

-- constructor
function Advection:new(tbl)
   local self = setmetatable({}, Advection)

   self._vel = ffi.new("double[7]", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) -- maximum 6D
   self._ndim = assert(tbl.ndim, "Eq.Advection: Must specify dimension using 'ndim'")
   assert(tbl.velocity, "Eq.Advection: Must specify velocity vector")

   -- read in velocity vector
   for d = 1, #tbl.velocity do
      self._vel[d] = tbl.velocity[d] 
   end

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Advection, { __call = function (self, o) return self.new(self, o) end })

-- Methods 
Advection.__index = {
   numEquations = function (self) return 1 end,
   numWaves = function (self) return 1 end,
   flux = function (self, dir, qIn, fOut)
      fOut[1] = self._vel[dir]*qIn[1]
   end,
   isPositive = function (self, q)
      return true
   end,
   rp = function (self, dir, delta, ql, qr, waves, s)
      waves[1][1] = delta[1]
      s[1] = self._vel[dir]
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
   Advection = Advection,
}


