-- Gkyl ------------------------------------------------------------------------
--
-- Burgers equation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"
local EqBase = require "Eq.EqBase"

-- Burgers equation
local Burgers = Proto(EqBase)

-- constructor
function Burgers:init(tbl)
   -- nothing to do
end

-- Methods 
function Burgers:numEquations() return 1 end
function Burgers:numWaves() return 1 end
function Burgers:flux(dir, qIn, fOut)
   fOut[1] = 0.5*qIn[1]*qIn[1]
end
function Burgers:isPositive(q)
   return true
end
function Burgers:rp(dir, delta, ql, qr, waves, s)
   waves[1][1] = delta[1]
   s[1] = 0.5*(ql[1]+qr[1])
end
function Burgers:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   amdq[1], apdq[1] = 0.0, 0.0
   if s[1] < 0 then
      amdq[1] = s[1]*waves[1][1]
   else
      apdq[1] = s[1]*waves[1][1]
   end
end

return Burgers


