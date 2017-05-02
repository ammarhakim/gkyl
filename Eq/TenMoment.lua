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

local _M = {}

-- C interfaces
ffi.cdef [[

/* TenMoment equations */
typedef struct {
   int64_t numFlux; /* Type of numerical flux to use */
   int useIntermediateWave; /* Flag to indicate intermediate wave use */
} TenMomentEqn_t;

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

-- Riemann problem for Ten-moment equations: `delta` is the vector we
-- wish to split, `ql`/`qr` the left/right states. On output, `waves`
-- and `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations. Note: This code is
-- essentially based on code used in my thesis i.e.
-- Miniwarpx. (A. Hakim)
local function rp(self, dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir`
end

-- The function to compute fluctuations is implemented as a template
-- which unrolls the inner loop
local qFluctuationsTempl = xsys.template([[
return function (dir, waves, s, amdq, apdq)
end
]])
-- function to compute fluctuations using q-wave method
local qFluctuations = loadstring( qFluctuationsTempl {} )()

local tenMoment_mt = {
   __new = function (self, tbl)
      local f = new(self)
      return f
   end,
   __index = {
      numEquations = function (self) return 10 end,
      numWaves = function (self) return 5 end,
      flux = function (self, dir, qIn, fOut)
	 local d = dirShuffle[dir] -- shuffle indices for `dir`
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
