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

-- helper to check if number if NaN
local function isNan(x) return x ~= x end

-- Riemann problem for Euler equations: `delta` is the vector we wish
-- to split, `ql`/`qr` the left/right states. On output, `waves` and
-- `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e. CLAWPACK and
-- Miniwarpx. (A. Hakim)
local function rp(self, dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir`
end

-- The function to compute fluctuations: fluctuations are hand-coded
-- from analytical expressions
local function qFluctuations(dir, waves, s, amdq, apdq)
end

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
