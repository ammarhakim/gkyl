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

local _M = {}

-- C interfaces
ffi.cdef [[

/* Euler equation for ideal gas */
typedef struct {
    double _gasGamma; /* Gas constant */
} EulerEqn_t;

void Euler_rp(EulerEqn_t *e, int dir, double *delta, double *ql, double *qr, double *waves, double *s);
]]

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int32_t[4]", 0, 2, 3, 4),
   new("int32_t[4]", 0, 3, 4, 2),
   new("int32_t[4]", 0, 4, 2, 3)
}

-- helper to check if number if NaN
local function isNan(x) return x ~= x end

-- Riemann problem for Euler equations: `delta` is the vector we wish
-- to split, `ql`/`qr` the left/right states. On output, `waves` and
-- `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e. CLAWPACK and
-- Miniwarpx. (A. Hakim)
local function rp(self, dir, delta, ql, qr, waves, s)
   ffi.C.Euler_rp(self, dir, delta:data(), ql:data(), qr:data(), waves:data(), s:data())
end

-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   local w1, w2, w3 = waves[1], waves[2], waves[3]
   local s1m, s2m, s3m = math.min(0, s[1]), math.min(0, s[2]), math.min(0, s[3])
   local s1p, s2p, s3p = math.max(0, s[1]), math.max(0, s[2]), math.max(0, s[3])

|for i = 1, 5 do
   amdq[${i}] = s1m*w1[${i}] + s2m*w2[${i}] + s3m*w3[${i}]
   apdq[${i}] = s1p*w1[${i}] + s2p*w2[${i}] + s3p*w3[${i}]

|end
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = loadstring( qFluctuationsTempl {} )()

local euler_mt = {
   __new = function (self, tbl)
      local f = new(self)
      f._gasGamma = assert(tbl.gasGamma, "Eq.Euler: Must specify gas adiabatic constant (gasGamma)")
      return f
   end,
   __index = {
      numEquations = function (self) return 5 end,
      numWaves = function (self) return 3 end,
      gasGamma = function (self) return self._gasGamma end,
      pressure = function (self, q)
	 return (self._gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
      end,
      flux = function (self, dir, qIn, fOut)
	 local d = dirShuffle[dir] -- shuffle indices for `dir`
	 local pr, u = self:pressure(qIn), qIn[d[1]]/qIn[1]
	 fOut[1] = qIn[d[1]] -- rho*u
	 fOut[d[1]] = qIn[d[1]]*u + pr -- rho*u*u + p
	 fOut[d[2]] = qIn[d[2]]*u -- rho*v*u
	 fOut[d[3]] = qIn[d[3]]*u -- rho*w*u
	 fOut[5] = (qIn[5]+pr)*u -- (E+p)*u
      end,
      isPositive = function (self, q)
	 if isNan(q[1]) or q[1] < 0.0 then return false end
	 local pr = self:pressure(q)
	 if isNan(pr) or pr < 0.0 then return false end
	 return true
      end,
      rp = rp,
      qFluctuations = function (self, dir, ql, qr, waves, s, amdq, apdq)
	 qFluctuations(dir, waves, s, amdq, apdq)
      end,
   }
}
_M.Euler = metatype(typeof("EulerEqn_t"), euler_mt)

return _M
