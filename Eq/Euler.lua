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
local Time = require "Lib.Time"

local _M = {}

-- C interfaces
ffi.cdef [[

/* Euler equation for ideal gas */
typedef struct {
    double _gasGamma; /* Gas constant */
   double _rpTime; /* Time spent in RP */
} EulerEqn_t;

]]

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int32_t[4]", 0, 2, 3, 4),
   new("int32_t[4]", 0, 3, 4, 2),
   new("int32_t[4]", 0, 4, 2, 3)
}

-- helper to check if number is NaN
local function isNan(x) return x ~= x end

-- Riemann problem for Euler equations: `delta` is the vector we wish
-- to split, `ql`/`qr` the left/right states. On output, `waves` and
-- `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e. CLAWPACK and
-- Miniwarpx. (A. Hakim)
local function rp(self, dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir`
   local g1 = self._gasGamma-1
   local rhol, rhor = ql[1], qr[1]
   local pl, pr = self:pressure(ql), self:pressure(qr)

   -- Roe averages: see Roe's original 1986 paper or LeVeque book
   local srrhol, srrhor = math.sqrt(rhol), math.sqrt(rhor)
   local ravgl1, ravgr1 = 1/srrhol, 1/srrhor
   local ravg2 = 1/(srrhol+srrhor)
   local u = (ql[d[1]]*ravgl1 + qr[d[1]]*ravgr1)*ravg2
   local v = (ql[d[2]]*ravgl1 + qr[d[2]]*ravgr1)*ravg2
   local w = (ql[d[3]]*ravgl1 + qr[d[3]]*ravgr1)*ravg2
   local enth = ((ql[5]+pl)*ravgl1 + (qr[5]+pr)*ravgr1)*ravg2   

   -- See http://ammar-hakim.org/sj/euler-eigensystem.html for
   -- notation and meaning of these terms
   local q2 = u*u+v*v+w*w
   local aa2 = g1*(enth-0.5*q2)
   local a = math.sqrt(aa2)
   local g1a2, euv = g1/aa2, enth-q2

   -- compute projections of jump
   local a4 = g1a2*(euv*delta[1] + u*delta[d[1]] + v*delta[d[2]] + w*delta[d[3]] - delta[5])
   local a2 = delta[d[2]] - v*delta[1]
   local a3 = delta[d[3]] - w*delta[1]
   local a5 = 0.5*(delta[d[1]] + (a-u)*delta[1] - a*a4)/a;
   local a1 = delta[1] - a4 - a5

   -- wave 1: eigenvalue is u-c
   local wv = waves[1]
   wv[1]  = a1
   wv[d[1]] = a1*(u-a)
   wv[d[2]] = a1*v
   wv[d[3]] = a1*w
   wv[5] = a1*(enth-u*a)
   s[1] = u-a

   -- wave 2: eigenvalue is u, u, u three waves are lumped into one
   wv = waves[2]
   wv[1]  = a4
   wv[d[1]] = a4*u
   wv[d[2]] = a4*v + a2
   wv[d[3]] = a4*w + a3
   wv[5] = a4*0.5*q2 + a2*v + a3*w
   s[2] = u

   -- wave 3: eigenvalue is u+c
   wv = waves[3]
   wv[1]  = a5
   wv[d[1]] = a5*(u+a)
   wv[d[2]] = a5*v
   wv[d[3]] = a5*w
   wv[5] = a5*(enth+u*a)
   s[3] = u+a
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
