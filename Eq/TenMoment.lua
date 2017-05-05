-- Gkyl ------------------------------------------------------------------------
--
-- Ten-moment equations in 3D
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

/* TenMoment equations */
typedef struct {
   int64_t numFlux; /* Type of numerical flux to use */
   int useIntermediateWave; /* Flag to indicate intermediate wave use */
   double _rpTime; /* Time spent in RP */
} TenMomentEqn_t;

]]

-- pre-defined constants to make life a little easier
local R, U, V, W, PXX, PXY, PXZ, PYY, PYZ, PZZ = 1, 1, 2, 3, 1, 2, 3, 4, 5, 6

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int32_t[4]", 0, 2, 3, 4),
   new("int32_t[4]", 0, 3, 4, 2),
   new("int32_t[4]", 0, 4, 2, 3)
}
local dirShufflePr = {
   new("int32_t[7]", 0, 5, 6, 7, 8, 9, 10),
   new("int32_t[7]", 0, 8, 9, 6, 10, 7, 5),
   new("int32_t[7]", 0, 10, 7, 9, 5, 6, 8)
}

-- helper to check if number is NaN
local function isNan(x) return x ~= x end

-- Multiply "w" by the phiPrime matrix (see
-- http://ammar-hakim.org/sj/tenmom-eigensystem.html). This brings a
-- right-eigenvector (or linear combinations) into conserved variable
-- space.
local function mulByPhiPrime(p0, u1, u2, u3, w, out) 
   out[1] = w[1] 
   out[2] = w[1]*u1+w[2]*p0 
   out[3] = w[1]*u2+w[3]*p0 
   out[4] = w[1]*u3+w[4]*p0 
   out[5] = w[1]*u1^2+2*w[2]*p0*u1+w[5] 
   out[6] = w[1]*u1*u2+w[2]*p0*u2+w[3]*p0*u1+w[6] 
   out[7] = w[1]*u1*u3+w[2]*p0*u3+w[4]*p0*u1+w[7] 
   out[8] = w[1]*u2^2+2*w[3]*p0*u2+w[8] 
   out[9] = w[1]*u2*u3+w[3]*p0*u3+w[4]*p0*u2+w[9] 
   out[10] = w[1]*u3^2+2*w[4]*p0*u3+w[10] 
end

-- function to convert conserved variables to primitive variables
local function primitive(q, out) 
   out[1] = q[1] 
   out[2] = q[2]/q[1] 
   out[3] = q[3]/q[1] 
   out[4] = q[4]/q[1] 
   out[5] = q[5]-q[2]^2/q[1] 
   out[6] = q[6]-(q[2]*q[3])/q[1] 
   out[7] = q[7]-(q[2]*q[4])/q[1] 
   out[8] = q[8]-q[3]^2/q[1] 
   out[9] = q[9]-(q[3]*q[4])/q[1] 
   out[10] = q[10]-q[4]^2/q[1] 
end

-- Riemann problem for Ten-moment equations: `delta` is the vector we
-- wish to split, `ql`/`qr` the left/right states. On output, `waves`
-- and `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my (Ammar) thesis i.e.
-- Miniwarpx. Parts of the code were generated automatically from
-- Maxima.
local function rp(self, dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir` for velocity
   local dp = dirShufflePr[dir] -- shuffle indices for `dir` for pressure tensor
   local sqrt = math.sqrt
   
   local vl, vr = ffi.new("double[11]"), ffi.new("double[11]")
   -- compute primitive values
   primitive(qr, vr)
   primitive(ql, vl)

   -- compute Roe averages
   local sqrl, sqrr = sqrt(vl[1]), sqrt(vr[1])
   local sqr1 = 1/(sqrl+sqrr)
      
   local p0 = sqrl*sqrr
   local u1 = (sqrl*vl[d[1]] + sqrr*vr[d[1]])*sqr1
   local u2 = (sqrl*vl[d[2]] + sqrr*vr[d[2]])*sqr1
   local u3 = (sqrl*vl[d[3]] + sqrr*vr[d[3]])*sqr1
   local p11 = (sqrr*vl[dp[1]]+sqrl*vr[dp[1]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[1]]-vl[d[1]])*(vr[d[1]]-vl[d[1]])
   local p12 = (sqrr*vl[dp[2]]+sqrl*vr[dp[2]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[1]]-vl[d[1]])*(vr[d[2]]-vl[d[2]])
   local p13 = (sqrr*vl[dp[3]]+sqrl*vr[dp[3]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[1]]-vl[d[1]])*(vr[d[3]]-vl[d[3]])
   local p22 = (sqrr*vl[dp[4]]+sqrl*vr[dp[4]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[2]]-vl[d[2]])*(vr[d[2]]-vl[d[2]])
   local p23 = (sqrr*vl[dp[5]]+sqrl*vr[dp[5]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[2]]-vl[d[2]])*(vr[d[3]]-vl[d[3]])
   local p33 = (sqrr*vl[dp[6]]+sqrl*vr[dp[6]])*sqr1 + 1.0/3.0*p0^2*sqr1^2*(vr[d[3]]-vl[d[3]])*(vr[d[3]]-vl[d[3]])

   -- local p0 = 0.5*(vl[1]+vr[1])
   -- local u1 = 0.5*(vl[d[1]]+vr[d[1]])
   -- local u2 = 0.5*(vl[d[2]]+vr[d[2]])
   -- local u3 = 0.5*(vl[d[3]]+vr[d[3]])
   -- local p11 = 0.5*(vl[dp[1]]+vr[dp[1]])
   -- local p12 = 0.5*(vl[dp[2]]+vr[dp[2]])
   -- local p13 = 0.5*(vl[dp[3]]+vr[dp[3]])
   -- local p22 = 0.5*(vl[dp[4]]+vr[dp[4]])
   -- local p23 = 0.5*(vl[dp[5]]+vr[dp[5]])
   -- local p33 = 0.5*(vl[dp[6]]+vr[dp[6]])

   local phiDelta = ffi.new("double[11]")   
   -- pre-multiply jump (delta) by phiPrime inverse: we do this as
   -- jumps are in conserved variables, while left eigenvectors used
   -- below are computed from primitive variables
   phiDelta[1] = delta[1] 
   phiDelta[2] = delta[d[1]]/p0-(delta[1]*u1)/p0 
   phiDelta[3] = delta[d[2]]/p0-(delta[1]*u2)/p0 
   phiDelta[4] = delta[d[3]]/p0-(delta[1]*u3)/p0 
   phiDelta[5] = delta[1]*u1^2-2*delta[d[1]]*u1+delta[dp[1]] 
   phiDelta[6] = delta[1]*u1*u2-delta[d[1]]*u2-delta[d[2]]*u1+delta[dp[2]] 
   phiDelta[7] = delta[1]*u1*u3-delta[d[1]]*u3-delta[d[3]]*u1+delta[dp[3]] 
   phiDelta[8] = delta[1]*u2^2-2*delta[d[2]]*u2+delta[dp[4]] 
   phiDelta[9] = delta[1]*u2*u3-delta[d[2]]*u3-delta[d[3]]*u2+delta[dp[5]] 
   phiDelta[10] = delta[1]*u3^2-2*delta[d[3]]*u3+delta[dp[6]]

   local leftProj = ffi.new("double[11]")
   -- project jumps on left eigenvectors (pray that LuaJIT eliminates common subexpressions)
   leftProj[1] = (phiDelta[2]*sqrt(p0)*p12)/(2*p11^(3/2))-(phiDelta[5]*p12)/(2*p11^2)-(phiDelta[3]*sqrt(p0))/(2*sqrt(p11))+phiDelta[6]/(2*p11) 
   leftProj[2] = (phiDelta[2]*sqrt(p0)*p13)/(2*p11^(3/2))-(phiDelta[5]*p13)/(2*p11^2)-(phiDelta[4]*sqrt(p0))/(2*sqrt(p11))+phiDelta[7]/(2*p11) 
   leftProj[3] = (-(phiDelta[2]*sqrt(p0)*p12)/(2*p11^(3/2)))-(phiDelta[5]*p12)/(2*p11^2)+(phiDelta[3]*sqrt(p0))/(2*sqrt(p11))+phiDelta[6]/(2*p11) 
   leftProj[4] = (-(phiDelta[2]*sqrt(p0)*p13)/(2*p11^(3/2)))-(phiDelta[5]*p13)/(2*p11^2)+(phiDelta[4]*sqrt(p0))/(2*sqrt(p11))+phiDelta[7]/(2*p11) 
   leftProj[5] = phiDelta[5]/(6*p11^2)-(phiDelta[2]*sqrt(p0))/(2*sqrt(3)*p11^(3/2)) 
   leftProj[6] = (phiDelta[2]*sqrt(p0))/(2*sqrt(3)*p11^(3/2))+phiDelta[5]/(6*p11^2) 
   leftProj[7] = phiDelta[1]-(phiDelta[5]*p0)/(3*p11) 
   leftProj[8] = (-(phiDelta[5]*p22)/(3*p11))+(4*phiDelta[5]*p12^2)/(3*p11^2)-(2*phiDelta[6]*p12)/p11+phiDelta[8] 
   leftProj[9] = (-(phiDelta[5]*p23)/(3*p11))+(4*phiDelta[5]*p12*p13)/(3*p11^2)-(phiDelta[6]*p13)/p11-(phiDelta[7]*p12)/p11+phiDelta[9] 
   leftProj[10] = (-(phiDelta[5]*p33)/(3*p11))+(4*phiDelta[5]*p13^2)/(3*p11^2)-(2*phiDelta[7]*p13)/p11+phiDelta[10]

   -- compute waves and speeds
   local wv = ffi.new("double[11]")

   -- Wave 1: (ev 1 and 2 are repeated)
   s[1] = u1-math.sqrt(p11/p0)
   wv[1] = 0.0 
   wv[d[1]] = 0.0 
   wv[d[2]] = -(1.0*leftProj[1]*sqrt(p11))/sqrt(p0) 
   wv[d[3]] = -(1.0*leftProj[2]*sqrt(p11))/sqrt(p0) 
   wv[dp[1]] = 0.0 
   wv[dp[2]] = leftProj[1]*p11 
   wv[dp[3]] = leftProj[2]*p11 
   wv[dp[4]] = 2.0*leftProj[1]*p12 
   wv[dp[5]] = leftProj[1]*p13+leftProj[2]*p12 
   wv[dp[6]] = 2.0*leftProj[2]*p13 

   mulByPhiPrime(p0, u1, u2, u3, wv, waves[1])
   
   -- Wave 2: (ev 3 and 4 are repeated)
   s[2] = u1+math.sqrt(p11/p0)
   wv[1] = 0.0 
   wv[d[1]] = 0.0 
   wv[d[2]] = (leftProj[3]*sqrt(p11))/sqrt(p0) 
   wv[d[3]] = (leftProj[4]*sqrt(p11))/sqrt(p0) 
   wv[dp[1]] = 0.0 
   wv[dp[2]] = leftProj[3]*p11 
   wv[dp[3]] = leftProj[4]*p11 
   wv[dp[4]] = 2.0*leftProj[3]*p12 
   wv[dp[5]] = leftProj[3]*p13+leftProj[4]*p12 
   wv[dp[6]] = 2.0*leftProj[4]*p13
   
   mulByPhiPrime(p0, u1, u2, u3, wv, waves[2])

   -- Wave 3 (ev 5)
   s[3] = u1-math.sqrt(3*p11/p0)
   wv[1] = leftProj[5]*p0*p11 
   wv[d[1]] = -(1.732050807568877*leftProj[5]*p11^(3/2))/sqrt(p0) 
   wv[d[2]] = -(1.732050807568877*leftProj[5]*sqrt(p11)*p12)/sqrt(p0) 
   wv[d[3]] = -(1.732050807568877*leftProj[5]*sqrt(p11)*p13)/sqrt(p0) 
   wv[dp[1]] = 3.0*leftProj[5]*p11^2 
   wv[dp[2]] = 3.0*leftProj[5]*p11*p12 
   wv[dp[3]] = 3.0*leftProj[5]*p11*p13 
   wv[dp[4]] = leftProj[5]*p11*p22+2.0*leftProj[5]*p12^2 
   wv[dp[5]] = leftProj[5]*p11*p23+2.0*leftProj[5]*p12*p13 
   wv[dp[6]] = leftProj[5]*p11*p33+2.0*leftProj[5]*p13^2
   
   mulByPhiPrime(p0, u1, u2, u3, wv, waves[3])
 
   -- Wave 4 (ev 6)
   s[4] = u1+math.sqrt(3*p11/p0)
   wv[1] = leftProj[6]*p0*p11 
   wv[d[1]] = (1.732050807568877*leftProj[6]*p11^(3/2))/sqrt(p0) 
   wv[d[2]] = (1.732050807568877*leftProj[6]*sqrt(p11)*p12)/sqrt(p0) 
   wv[d[3]] = (1.732050807568877*leftProj[6]*sqrt(p11)*p13)/sqrt(p0) 
   wv[dp[1]] = 3.0*leftProj[6]*p11^2 
   wv[dp[2]] = 3.0*leftProj[6]*p11*p12 
   wv[dp[3]] = 3.0*leftProj[6]*p11*p13 
   wv[dp[4]] = leftProj[6]*p11*p22+2.0*leftProj[6]*p12^2 
   wv[dp[5]] = leftProj[6]*p11*p23+2.0*leftProj[6]*p12*p13 
   wv[dp[6]] = leftProj[6]*p11*p33+2.0*leftProj[6]*p13^2 

   mulByPhiPrime(p0, u1, u2, u3, wv, waves[4])

   -- Wave 5: (ev 7, 8, 9, 10 are repeated)
   s[5] = u1
   wv[1] = leftProj[7] 
   wv[d[1]] = 0 
   wv[d[2]] = 0 
   wv[d[3]] = 0 
   wv[dp[1]] = 0 
   wv[dp[2]] = 0 
   wv[dp[3]] = 0 
   wv[dp[4]] = leftProj[8] 
   wv[dp[5]] = leftProj[9] 
   wv[dp[6]] = leftProj[10]    

   mulByPhiPrime(p0, u1, u2, u3, wv, waves[5])
end

-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   local w1, w2, w3, w4, w5 = waves[1], waves[2], waves[3], waves[4], waves[5]
   local s1m, s2m, s3m, s4m, s5m = math.min(0, s[1]), math.min(0, s[2]), math.min(0, s[3]), math.min(0, s[4]), math.min(0, s[5])
   local s1p, s2p, s3p, s4p, s5p = math.max(0, s[1]), math.max(0, s[2]), math.max(0, s[3]), math.max(0, s[4]), math.max(0, s[5])

|for i = 1, 10 do
   amdq[${i}] = s1m*w1[${i}] + s2m*w2[${i}] + s3m*w3[${i}] + s4m*w4[${i}] + s5m*w5[${i}]
   apdq[${i}] = s1p*w1[${i}] + s2p*w2[${i}] + s3p*w3[${i}] + s4p*w4[${i}] + s5p*w5[${i}]

|end
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = loadstring( qFluctuationsTempl {} )()

local tenMoment_mt = {
   __new = function (self, tbl)
      local f = new(self)
      f._rpTime = 0.0
      return f
   end,
   __index = {
      numEquations = function (self) return 10 end,
      numWaves = function (self) return 5 end,
      flux = function (self, dir, qIn, fOut)
	 local d = dirShuffle[dir] -- shuffle indices for `dir`
	 local dp = dirShufflePr[dir] -- shuffle indices for `dir` for pressure tensor
	 local v = ffi.new("double[11]")
	 
	 primitive(qIn, v)
	 
	 fOut[1] = qIn[d[1]]
	 fOut[d[1]] = qIn[dp[PXX]]
	 fOut[d[2]] = qIn[dp[PXY]]
	 fOut[d[3]] = qIn[dp[PXZ]]
	 fOut[dp[1]] = v[R]*v[d[U]]^3 + 3*v[d[U]]*v[dp[PXX]]
	 fOut[dp[2]] = v[R]*v[d[U]]^2*v[d[V]] + 2*v[d[U]]*v[dp[PXY]] + v[d[V]]*v[dp[PXX]]
	 fOut[dp[3]] = v[R]*v[d[U]]^2*v[d[W]] + 2*v[d[U]]*v[dp[PXZ]] + v[d[W]]*v[dp[PXX]]
	 fOut[dp[4]] = v[R]*v[d[U]]*v[d[V]]^2 + v[d[U]]*v[dp[PYY]] + 2*v[d[V]]*v[dp[PXY]]
	 fOut[dp[5]] = v[R]*v[d[U]]*v[d[V]]*v[d[W]] + v[d[U]]*v[dp[PYZ]] + v[d[V]]*v[dp[PXZ]] + v[d[W]]*v[dp[PXY]]
	 fOut[dp[6]] = v[R]*v[d[U]]*v[d[W]]^2 + v[d[U]]*v[dp[PZZ]] + 2*v[d[W]]*v[dp[PXZ]]
      end,
      isPositive = function (self, q)
	 return true
      end,
      rp = rp,
      qFluctuations = function (self, dir, ql, qr, waves, s, amdq, apdq)
	 qFluctuations(dir, waves, s, amdq, apdq)
      end,
   }
}
_M.TenMoment = metatype(typeof("TenMomentEqn_t"), tenMoment_mt)

return _M
