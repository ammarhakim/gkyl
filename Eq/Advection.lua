-- Gkyl ------------------------------------------------------------------------
--
-- Linear, constant velocity advection equation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local Proto = require "Lib.Proto"
local EqBase = require "Eq.EqBase"

-- Constant velocity advection equation in upto NDIM=6
local Advection = Proto(EqBase)

-- constructor
function Advection:init(tbl)
   self._vel = ffi.new("double[7]", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) -- maximum 6D
   assert(tbl.velocity, "Eq.Advection: Must specify velocity vector")
   self._ndim = #tbl.velocity

   -- read in velocity vector
   for d = 1, #tbl.velocity do
      self._vel[d] = tbl.velocity[d] 
   end
end

-- Methods 
function Advection:numEquations() return 1 end
function Advection:numWaves() return 1 end
function Advection:velocity(dir) return self._vel[dir] end

function Advection:flux(dir, qIn, fOut)
   fOut[1] = self._vel[dir]*qIn[1]
end

function Advection:isPositive(q) return true end

function Advection:rp(dir, delta, ql, qr, waves, s)
   waves[1][1] = delta[1]
   s[1] = self._vel[dir]
end

function Advection:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   amdq[1], apdq[1] = 0.0, 0.0
   if s[1] < 0 then
      amdq[1] = s[1]*waves[1][1]
   else
      apdq[1] = s[1]*waves[1][1]
   end
end

return Advection
