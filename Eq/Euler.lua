-- Gkyl ------------------------------------------------------------------------
--
-- Euler (ideal gas) equations in 3D
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

_M = {}

-- C interfaces
ffi.cdef [[

/* Euler equation for ideal gas */
typedef struct {
    double _gasGamma; /* Gas constant */
} EulerEqn_t;

]]

-- helper to check if number if NaN
local function isNan(x) return x ~= x end

local euler_mt = {
   __new = function (self, tbl)
      local f = new(self)
      f._gasGamma = assert(tbl.gasGamma, "Eq.Euler: Must specify gas adiabatic constant (gasGamma)")
      return f
   end,
   __index = {
      numEquations = function (self) return 5 end,
      numWaves = function (self) return 3 end,
      gasGamma = function (self)
	 return self._gasGamma
      end,
      pressure = function (self, q)
	 return (self._gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
      end,
      flux = function (self, qIn, fOut)
	 local rho1, pr = 1/qIn[1], self:pressure(qIn)
	 local u = qIn[2]*rho1

	 fOut[1] = qIn[2] -- rho*u
	 fOut[2] = qIn[2]*u + pr -- rho*u*u + p
	 fOut[3] = qIn[3]*u -- rho*v*u
	 fOut[4] = qIn[4]*u -- rho*w*u
	 fOut[5] = (qIn[5]+pr)*u -- (E+p)*u
      end,
      isPositive = function (self, q)
	 if isNan(q[1]) or q[1] < 0.0 then return false end
	 local pr = self:pressure(q)
	 if isNan(pr) or pr < 0.0 then return false end
	 return true
      end,
   }
}
_M.Euler = metatype(typeof("EulerEqn_t"), euler_mt)

return _M
