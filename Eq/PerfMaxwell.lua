-- Gkyl ------------------------------------------------------------------------
--
-- Perfectly Hyperbolic Maxwell equations. See
-- http://ammar-hakim.org/sj/maxwell-eigensystem.html and references
-- on that page for details.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local BoundaryCondition = require "Updater.BoundaryCondition"
local EqBase = require "Eq.EqBase"
local MaxwellModDecl = require "Eq.maxwellData.MaxwellModDecl"
local Proto = require "Lib.Proto"
local ffi = require "ffi"
local xsys = require "xsys"

-- Resuffle indices for various direction Riemann problem. The first
-- entry is just a buffer to allow 1-based indexing
local dirShuffle = {
   ffi.new("int32_t[7]", 0, 1, 2, 3, 4, 5, 6),
   ffi.new("int32_t[7]", 0, 2, 3, 1, 5, 6, 4),
   ffi.new("int32_t[7]", 0, 3, 1, 2, 6, 4, 5)
}

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

-- Perfectly hyperbolic Maxwell equations
local PerfMaxwell = Proto(EqBase)

-- ctor
function PerfMaxwell:init(tbl)
   self._c = assert(tbl.lightSpeed, "Eq.PerfMaxwell: Must specify light speed (lightSpeed)")
   self._ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0

   -- if specified, store basis functions
   self._basis = tbl.basis and tbl.basis or nil

   -- tau parameter used for adding extra (less) diffusion to Ampere-Maxwell, while adding less (more) diffusion to Faraday equation
   -- if no tau parameter is specified, defaults to the speed of light
   self._tau = tbl.tau and tbl.tau or self._c

   -- store pointers to C kernels implementing volume and surface
   -- terms
   self._volTerm, self._surfTerms = nil, nil
   if self._basis then
      local nm, ndim, p = self._basis:id(), self._basis:ndim(), self._basis:polyOrder()
      self._volTerm = MaxwellModDecl.selectVol(nm, ndim, p)

      -- numFlux used for selecting which type of numerical flux function to use
      -- default is "upwind," supported options: "central," "upwind"
      self._numFlux = tbl.numFlux and tbl.numFlux or "upwind"
      if self._numFlux == "upwind" then
         self._surfTerms = MaxwellModDecl.selectSurf(nm, ndim, p)
      elseif self._numFlux == "central" then
         print("selecting central fluxes")
         self._surfTerms = MaxwellModDecl.selectCentralSurf(nm, ndim, p)
      else
         assert(self._numFLux, "Eq.PerfMaxwell: Incorrect numerical flux specified, options supported: 'central' and 'upwind' ")
      end
   end

   -- maximum characteristic speed
   self._maxc = math.max(self._c, self._c*self._ce, self._c*self._cb)

   -- store stuff in C struct for use in DG solvers
   self._ceqn = ffi.new("MaxwellEq_t", {self._c, self._ce, self._cb})
end

-- Methods
function PerfMaxwell:numEquations() return 8 end
function PerfMaxwell:numWaves() return 6 end
function PerfMaxwell:isPositive(q) return true end

-- flux in direction dir
function PerfMaxwell:flux(dir, qIn, fOut)
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
end

-- Riemann problem for Perfectly Hyperbolic Maxwell equations: `delta`
-- is the vector we wish to split, `ql`/`qr` the left/right states. On
-- output, `waves` and `s` contain the waves and speeds. waves is a
-- mwave X meqn matrix. See LeVeque's book for explanations. Note:
-- This code is essentially based on code used in my thesis
-- i.e. CLAWPACK and Miniwarpx. (A. Hakim)
function PerfMaxwell:rp(dir, delta, ql, qr, waves, s)
   local d = dirShuffle[dir] -- shuffle indices for `dir`
   local c, c1 = self._c, 1/self._c
   local v = c
   local a = 0.5 * (v/c - 1)

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
   ffi.fill(waves:data(), 8*6*ffi.sizeof("double"))

   -- wave 1:
   local w = waves[1]
   w[d[4]] = a1
   w[8] = -a1*c
   s[1] = -c*self._cb

   -- wave 2:
   local w = waves[2]
   w[d[4]] = a2
   w[8] = a2*c
   s[2] = c*self._cb

   -- wave 3:
   local w = waves[3]
   w[d[1]] = a3
   w[7] = -a3*c1
   s[3] = -c*self._ce

   -- wave 4:
   local w = waves[4]
   w[d[1]] = a4
   w[7] = a4*c1
   s[4] = c*self._ce

   -- wave 5: (two waves with EV -c, -c lumped into one)
   local w = waves[5]
   w[d[2]] = a5
   w[d[3]] = a6
   w[d[5]] = a6*c1 + a * delta[d[5]]
   w[d[6]] = -a5*c1 + a * delta[d[6]]
   s[5] = -c

   -- wave 6: (two waves with EV c, c lumped into one)
   local w = waves[6]
   w[d[2]] = a7
   w[d[3]] = a8
   w[d[5]] = -a8*c1 + a * delta[d[5]]
   w[d[6]] = a7*c1 + a * delta[d[6]]
   s[6] = c
end

-- Compute q-fluctuations
function PerfMaxwell:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   qFluctuations(dir, waves, s, amdq, apdq)

   local c = s[6]
   local v = c
   local d = dirShuffle[dir] -- shuffle indices for `dir`
   -- for _,i in ipairs({2, 3, 5, 6}) do
   for _,i in ipairs({5, 6}) do
      apdq[d[i]] = apdq[d[i]] + 0.5 * (v-c) * (qr[d[i]] - ql[d[i]])
      amdq[d[i]] = amdq[d[i]] - 0.5 * (v-c) * (qr[d[i]] - ql[d[i]])
   end
   local v = c
   local d = dirShuffle[dir] -- shuffle indices for `dir`
   -- for _,i in ipairs({2, 3, 5, 6}) do
   -- for _,i in ipairs({5, 6}) do
   for _,i in ipairs({2, 3}) do
      apdq[d[i]] = apdq[d[i]] + 0.5 * (v-c) * (qr[d[i]] - ql[d[i]])
      amdq[d[i]] = amdq[d[i]] - 0.5 * (v-c) * (qr[d[i]] - ql[d[i]])
   end

end

-- Maximum wave speed
function PerfMaxwell:maxSpeed(dir, w, dx, q)
   return self._maxc
end       
-- Volume integral term for use in DG scheme
function PerfMaxwell:volTerm(w, dx, idx, q, out)
   return self._volTerm(self._ceqn, w:data(), dx:data(), q:data(), out:data())
end
-- Surface integral term for use in DG scheme
function PerfMaxwell:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   return self._surfTerms[dir](self._ceqn, wl:data(), wr:data(), dxl:data(), dxr:data(), self._tau, ql:data(), qr:data(), outl:data(), outr:data())
end

-- Create and add BCs specific to Maxwell equations

-- PEC
local bcCondWallElc = BoundaryCondition.ZeroTangent { components = {1, 2, 3} }
local bcCondWallMgn = BoundaryCondition.ZeroNormal { components = {4, 5, 6} }
local bcCondWallPot = BoundaryCondition.Copy { components = {7, 8}, fact = {-1, 1} }
PerfMaxwell.bcCondWall = { bcCondWallElc, bcCondWallMgn, bcCondWallPot  }

function PerfMaxwell:setAuxFields(auxFields)
end

return PerfMaxwell
