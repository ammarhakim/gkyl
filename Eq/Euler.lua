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

-- functions to rotate conserved variables to local aligned
-- coordinates
alignedRotToLocal = {}
alignedRotToLocal[1] = function (qin, qout)
   qout[1], qout[2], qout[3], qout[4], qout[5] = qin[1], qin[2], qin[3], qin[4], qin[5]
end
alignedRotToLocal[2] = function (qin, qout)
   qout[1], qout[2], qout[3], qout[4], qout[5] = qin[1], qin[3], -qin[2], qin[4], qin[5]
end
alignedRotToLocal[3] = function(qin, qout)
   error("Not implemented")
end

-- functions to rotate conserved variables to global aligned
-- coordinates
alignedRotToGlobal = {}
alignedRotToGlobal[1] = function (qin, qout)
   qout[1], qout[2], qout[3], qout[4], qout[5] = qin[1], qin[2], qin[3], qin[4], qin[5]
end
alignedRotToGlobal[2] = function (qin, qout)
   qout[1], qout[2], qout[3], qout[4], qout[5] = qin[1], -qin[3], qin[2], qin[4], qin[5]
end
alignedRotToGlobal[3] = function (qin, qout)
   error("Not implemented")
end

-- Riemann problem for Euler equations: `delta` is the vector we wish
-- to split, `ql`/`qr` the left/right states. On output, `waves` and
-- `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e. CLAWPACK and
-- Miniwarpx. (A. Hakim)
local function rp(self, delta, ql, qr, waves, s)
   local g1 = self._gasGamma-1
   local rhol, rhor = ql[1], qr[1]
   local pl, pr = self:pressure(ql), self:pressure(qr)

   -- Roe averages: see Roe's original 1986 paper or LeVeque book
   local ravgl1, ravgr1 = 1/math.sqrt(rhol), 1/math.sqrt(rhor)
   local ravg2 = 1/(math.sqrt(rhol)+math.sqrt(rhor))
   local u = (ql[2]*ravgl1 + qr[2]*ravgr1)*ravg2
   local v = (ql[3]*ravgl1 + qr[3]*ravgr1)*ravg2
   local w = (ql[4]*ravgl1 + qr[4]*ravgr1)*ravg2
   local enth = ((ql[5]+pl)*ravgl1 + (qr[5]+pr)*ravgr1)*ravg2   

    -- See http://ammar-hakim.org/sj/euler-eigensystem.html for
    -- notation and meaning of these terms
   local q2 = u*u+v*v+w*w
   local aa2 = g1*(enth-0.5*q2)
   local a = math.sqrt(aa2)
   local g1a2, euv = g1/aa2, enth-q2

   -- compute projections of jump
   local a4 = g1a2*(euv*delta[1] + u*delta[2] + v*delta[3] + w*delta[4] - delta[5])
   local a2 = delta[3] - v*delta[1]
   local a3 = delta[4] - w*delta[1]
   local a5 = 0.5*(delta[2] + (a-u)*delta[1] - a*a4)/a;
   local a1 = delta[1] - a4 - a5

   -- wave 1: eigenvalue is u-c
   local wv = waves[1]
   wv[1]  = a1
   wv[2] = a1*(u-a)
   wv[3] = a1*v
   wv[4] = a1*w
   wv[5] = a1*(enth-u*a)
   s[1] = u-a

   -- wave 2: eigenvalue is u, u, u three waves are lumped into one
   wv = waves[2]
   wv[1]  = a4
   wv[2] = a4*u
   wv[3] = a4*v + a2
   wv[4] = a4*w + a3
   wv[5] = a4*0.5*q2 + a2*v + a3*w
   s[2] = u

   -- wave 3: eigenvalue is u+c
   wv = waves[3]
   wv[1]  = a5
   wv[2] = a5*(u+a)
   wv[3] = a5*v
   wv[4] = a5*w
   wv[5] = a5*(enth+u*a)
   s[3] = u+a
end

-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (waves, s, amdq, apdq)
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
      flux = function (self, qIn, fOut)
	 local pr, u = self:pressure(qIn), qIn[2]/qIn[1]
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
      rp = rp,
      qFluctuations = function (self, ql, qr, waves, s, amdq, apdq)
	 return qFluctuations(waves, s, amdq, apdq)
      end,
      rotateToLocalAligned = function (self, dir, qin, qout)
	 alignedRotToLocal[dir](qin, qout)
      end,
      rotateToGlobalAligned = function (self, dir, qin, qout)
	 alignedRotToGlobal[dir](qin, qout)
      end,
   }
}
_M.Euler = metatype(typeof("EulerEqn_t"), euler_mt)

return _M
