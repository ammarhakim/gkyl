-- Gkyl ------------------------------------------------------------------------
--
-- Ten-moment equations in 3D.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")
local Time = require "Lib.Time"
local BoundaryCondition = require "Updater.BoundaryCondition"

local _M = {}

-- C interfaces.
ffi.cdef [[

/* TenMoment equations */
typedef struct {
   int _numWaves; /* Number of waves in system */
   int _rpType; /* Type of RP to use */
   int64_t numFlux; /* Type of numerical flux to use */
   int useIntermediateWave; /* Flag to indicate intermediate wave use */
   double _rpTime; /* Time spent in RP */
   double _tv[11]; /* To store stuff in flux() function */
} TenMomentEqn_t;

/* Ten-moment RP solver */
  void gkylTenMomentRp(int dir, double delta[10], double ql[10], double qr[10], double *waves, double s[5]);
/* Ten-moment q-fluctuations */
  void gkylTenMomentQFluct(const int dir, double *ql, double *qr, double *waves, double *s, double *amdq, double *apdq);

/* Ten-moment RP solver using Lax flux */
  void gkylTenMomentRpLax(int dir, double delta[10], double ql[10], double qr[10], double *waves, double s[5]);
/* Ten-moment q-fluctuations for Lax flux */
  void gkylTenMomentQFluctLax(const int dir, double *ql, double *qr, double *waves, double *s, double *amdq, double *apdq);
]]

-- Pre-defined constants to make life a little easier.
local R, U, V, W, PXX, PXY, PXZ, PYY, PYZ, PZZ = 1, 1, 2, 3, 1, 2, 3, 4, 5, 6

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   new("int[4]", 0, 2, 3, 4),
   new("int[4]", 0, 3, 4, 2),
   new("int[4]", 0, 4, 2, 3)
}
local dirShufflePr = {
   new("int[7]", 0, 5, 6, 7, 8, 9, 10),
   new("int[7]", 0, 8, 9, 6, 10, 7, 5),
   new("int[7]", 0, 10, 7, 9, 5, 6, 8)
}

-- helper to check if number is NaN
local function isNan(x) return x ~= x end

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

local RP_TYPE_ROE, RP_TYPE_LAX  = 1, 2

-- Riemann problem for Ten-moment equations: `delta` is the vector we
-- wish to split, `ql`/`qr` the left/right states. On output, `waves`
-- and `s` contain the waves and speeds. waves is a mwave X meqn
-- matrix. See LeVeque's book for explanations.
local function rp(self, dir, delta, ql, qr, waves, s)
   -- dispatch call to C++ function implementing the RP
   if self._rpType == RP_TYPE_ROE then
      ffi.C.gkylTenMomentRp(dir, delta:data(), ql:data(), qr:data(), waves:data(), s:data())
   else
      ffi.C.gkylTenMomentRpLax(dir, delta:data(), ql:data(), qr:data(), waves:data(), s:data())
   end
end

local tenMoment_mt = {
   __new = function (self, tbl)
      local f = new(self)
      f._rpType = RP_TYPE_ROE
      if tbl.numericalFlux then
	 if tbl.numericalFlux == "roe" then
	    f._rpType = RP_TYPE_ROE
	 elseif tbl.numericalFlux == "lax" then
	    f._rpType = RP_TYPE_LAX
	 else
	    assert(false, "Euler: 'numericalFlux' must be one of 'roe' or 'lax'")
	 end
      end
      f._numWaves = f._rpType == RP_TYPE_ROE and 5 or 1 -- Number of waves.
      f._rpTime = 0.0
      return f
   end,
   __index = {
      numEquations = function (self) return 10 end,
      numWaves = function (self) return self._numWaves end,
      flux = function (self, dir, qIn, fOut)
	 local d = dirShuffle[dir] -- shuffle indices for `dir`
	 local dp = dirShufflePr[dir] -- shuffle indices for `dir` for pressure tensor
	 local v = self._tv
	 
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
          if isNan(q[1]) or q[1] < 0.0 then return false end
          local pxx = q[5] - q[2]^2 / q[1]
          if isNan(pxx) or pxx < 0.0 then return false end
          local pyy = q[8] - q[3]^2 / q[1]
          if isNan(pyy) or pyy < 0.0 then return false end
          local pzz = q[10] - q[4]^2 / q[1]
          if isNan(pzz) or pzz < 0.0 then return false end
          return true
      end,
      rp = rp,
      qFluctuations = function (self, dir, ql, qr, waves, s, amdq, apdq)
	 if self._rpType == RP_TYPE_ROE then
	    ffi.C.gkylTenMomentQFluct(dir, ql:data(), qr:data(), waves:data(),
                                      s:data(), amdq:data(), apdq:data())
         else
	    ffi.C.gkylTenMomentQFluctLax(dir, ql:data(), qr:data(),
                                         waves:data(), s:data(), amdq:data(),
                                         apdq:data())
         end
      end,
   }
}
local TenMomentObj = metatype(typeof("TenMomentEqn_t"), tenMoment_mt)

-- create a wrapper on Ten-moment eqn object and provide BCs specific
-- to equations
local TenMoment = {}
function TenMoment:new(tbl)
   local self = setmetatable({}, TenMoment)
   return TenMomentObj(tbl)
end
-- make object callable, and redirect call to the :new method
setmetatable(TenMoment, { __call = function (self, o) return self.new(self, o) end })

-- rho, Pxx, Pyy, Pzz
local bcWallCopy = BoundaryCondition.Copy { components = {1, 5, 8, 10} }
local bcWallFlow = BoundaryCondition.ZeroNormal { components = {2, 3, 4} }
-- Pyz, Pxz, Pxy
local bcWallPres = BoundaryCondition.ZeroTangent { components = {9, 7, 6} }
TenMoment.bcWall = { bcWallCopy, bcWallFlow, bcWallPres }

TenMoment.bcConst = function(rho, rhou, rhov, rhow, pxx, pxy, pxz, pyy, pyz, pzz)
   local bc = BoundaryCondition.Const {
      components = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
      values     = {rho, rhou, rhov, rhow, pxx, pxy, pxz, pyy, pyz, pzz}
   }
   return { bc }
end

return TenMoment
