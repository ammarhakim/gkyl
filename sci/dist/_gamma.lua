--------------------------------------------------------------------------------
-- Gamma statistical distribution.
--
-- Copyright (C) 2011-2015 Stefano Peluchetti. All rights reserved.
--------------------------------------------------------------------------------

-- TODO: consider branchless-free loop (see Mike mail) for vector sampling, 
-- TODO: may be used in any case via an (heuristic-adaptive) buffer-based algo.

-- TODO: consider using a different metatable for the two algorithms.

local xsys    = require "xsys"
local ffi     = require "ffi"
local math    = require "sci.math"
local _normal = require "sci.dist._normal"

local exp, log, sqrt, min, gamma, loggamma = xsys.from(math, 
     "exp, log, sqrt, min, gamma, loggamma")

local nd01 = _normal.dist(0, 1)

-- Based on Marsaglia-Tsang method, valid for a >= 1, see:
-- http://www.cparity.com/projects/AcmClassification/samples/358414.pdf .
-- PERF: not worth checking for v <= 0 separately.
local function samplea1b(a, b, rng)
  local d = a - 1/3
  local c = 1/sqrt(9*d)
  
  do -- PERF: avoid (>95% of the time) problematic loop below with this: 
    local x = nd01:sample(rng)
    local v = (1 + c*x)^3
    local u = rng:sample()
    if min(v, 1 - 0.0331*x^4 - u) > 0 -- Squeeze, around 92% of the time.
    or min(v, -log(u) + 0.5*x^2 + d*(1 - v + log(v))) > 0 then
      return v*d/b -- v*d ~ gamma(a), scaling: v*d/b ~ gamma(a, b).
    end
  end
  
  -- PERF: better this than 'repeat until' alternative, as it allows to return
  -- PERF: directly from within the loop:
  while true do
    local x = nd01:sample(rng)
    local v = (1 + c*x)^3
    -- PERF: squeeze step not worth here.
    if min(v, -log(rng:sample()) + 0.5*x^2 + d*(1 - v + log(v))) > 0 then
      return v*d/b -- v*d ~ gamma(a), scaling: v*d/b ~ gamma(a, b).
    end
  end
end

local function sampleab(a, b, rng)
  if a >= 1 then
    return samplea1b(a, b, rng)
  else
    local y = samplea1b(a + 1, b, rng)
    return rng:sample()^(1/a)*y
  end
end
     
local gamm_mt = {
  __new = function(ct, alpha, beta)
    if not alpha or not beta then
      error("distribution parameters must be set at construction")
    end
    if alpha <= 0 or beta <= 0 then
      error("alpha and beta must be positive, alpha is "..alpha..", beta is "
        ..beta)
    end
    return ffi.new(ct, alpha, beta)
  end,
  copy = function(self)
    return ffi.new(ffi.typeof(self), self)
  end,
  range = function(self)
    return 0, 1/0
  end,
  pdf = function(self, x)
    if x < 0 then return 0 end
    local a, b = self._a, self._b
    return b^a * x^(a - 1) * exp(-b*x) / gamma(a)
  end,
  logpdf = function(self, x)
    if x < 0 then return -1/0 end
    local a, b = self._a, self._b
    return a*log(b) - loggamma(a) + (a - 1)*log(x) - b*x
  end,
  mean = function(self)
    return self._a/self._b
  end,
  var = function(self)
    return self._a/self._b^2
  end,
  sample = function(self, rng)
    return sampleab(self._a, self._b, rng)
  end,
}
gamm_mt.__index = gamm_mt

local dist = ffi.metatype("struct { double _a, _b; }", gamm_mt)

return {
  dist     = dist,
  sampleab = sampleab,
}