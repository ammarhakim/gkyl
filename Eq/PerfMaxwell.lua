-- Gkyl ------------------------------------------------------------------------
--
-- Perfectly Hyperbolic Maxwell equations. See
-- http://ammar-hakim.org/sj/maxwell-eigensystem.html and references
-- on that page for details.
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

/* PH Maxwell equations */
typedef struct {
    double _c; /* Speed of light */
    double _ce, _cb; /* Electric and Magnetic field error speeds (multiple of _c) */
} PerfMaxwellEqn_t;

]]

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int32_t[7]", 0, 1, 2, 3, 4, 5, 6),
   new("int32_t[7]", 0, 2, 3, 1, 5, 6, 4),
   new("int32_t[7]", 0, 3, 1, 2, 6, 4, 5)
}

-- Riemann problem for Euler equations: `delta` is the vector we wish
-- to split, `ql`/`qr` the left/right states. On output, `waves` and
-- `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e. CLAWPACK and
-- Miniwarpx. (A. Hakim)
local function rp(self, dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir`
   local c, c1 = self._c, 1/self._c

   -- compute projections of jump (generated from Maxima)
   local a1 = 0.5*(delta[d[4]]-delta[8]*c1)
   local a2 = 0.5*(delta[8]*c1+delta[d[4]])
   local a3 = 0.5*(delta[d[1]]-delta[7]*c)
   local a4 = 0.5*(delta[7]*c+delta[d[1]])
   local a5 = 0.5*(delta[d[2]]-delta[d[6]]*c)
   local a6 = 0.5*(delta[d[5]]*c+delta[d[3]])
   local a7 = 0.5*(delta[d[6]]*c+delta[d[2]])
   local a8 = 0.5*(delta[d[3]]-delta[d[5]]*c)

   -- set waves to 0.0 as most entries vanish
   fill(waves:data(), 8*6*sizeof("double"))

   -- wave 1:
   waves[1][d[4]] = a1
   waves[1][8] = -a1*c
   s[1] = -c*self._cb

   -- wave 2:
   waves[2][d[4]] = a2
   waves[2][8] = a2*c
   s[2] = c*self._cb

   -- wave 3:
   waves[3][d[1]] = a3
   waves[3][7] = -a3*c1
   s[3] = -c*self._ce

   -- wave 4:
   waves[4][d[1]] = a4
   waves[4][7] = a4*c1
   s[4] = c*self._ce

   -- wave 5: (two waves with EV -c, -c lumped into one)
   waves[5][d[2]] = a5
   waves[5][d[3]] = a6
   waves[5][d[5]] = a6*c1
   waves[5][d[6]] = -a5*c1
   s[5] = -c

   -- wave 6: (two waves with EV c, c lumped into one)
   waves[6][d[2]] = a7
   waves[6][d[3]] = a8
   waves[6][d[5]] = -a8*c1
   waves[6][d[6]] = a7*c1
   s[6] = c
end

-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
   local w1p, w2p, w3p = waves[2], waves[4], waves[6]
   local s1p, s2p, s3p = s[2], s[4], s[6]

|for i = 1, 8 do
   apdq[${i}] = s1p*w1p[${i}] + s2p*w2p[${i}] + s3p*w3p[${i}]
|end

   local w1m, w2m, w3m = waves[1], waves[3], waves[5]
   local s1m, s2m, s3m = s[1], s[3], s[5]

|for i = 1, 8 do
   amdq[${i}] = s1m*w1m[${i}] + s2m*w2m[${i}] + s3m*w3m[${i}]
|end
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = loadstring( qFluctuationsTempl {} )()

local maxwell_mt = {
   __new = function (self, tbl)
      local f = new(self)
      f._c = assert(tbl.lightSpeed, "Eq.PerfMaxwell: Must specify gas light speed (lightSpeed)")
      f._ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
      f._cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0
      return f
   end,
   __index = {
      numEquations = function (self) return 8 end,
      numWaves = function (self) return 6 end,
      flux = function (self, dir, qIn, fOut)
	 local d = dirShuffle[dir] -- shuffle indices for `dir`
	 local c2 = self._c*self._c

	 fOut[d[1]] = self._ce*c2*qIn[7]
	 fOut[d[2]] = c2*qIn[d[6]]
	 fOut[d[3]] = -c2*qIn[d[5]]
	 fOut[d[4]] = self._cb*qIn[8]
	 fOut[d[5]] = -qIn[d[3]]
	 fOut[d[6]] = qIn[d[2]]
	 fOut[7] = self._ce*qIn[d[1]]
	 fOut[8] = self._cb*c2*qIn[d[4]]
      end,
      isPositive = function (self, q)
	 return true
      end,
      rp = rp,
      qFluctuations = function (self, dir, ql, qr, waves, s, amdq, apdq)
	 return qFluctuations(dir, waves, s, amdq, apdq)
      end,
   }
}
_M.PhMaxwell = metatype(typeof("PerfMaxwellEqn_t"), maxwell_mt)

return _M
