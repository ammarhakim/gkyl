-- Gkyl ------------------------------------------------------------------------
--
-- App support code: VlasovProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

--local Time = require "Lib.Time".
local FunctionProjectionParent   = require ("App.Projection.KineticProjection").FunctionProjection
local MaxwellianProjectionParent = require ("App.Projection.KineticProjection").MaxwellianProjection
local Proto                      = require "Lib.Proto"

----------------------------------------------------------------------
-- Vlasov-specific VlasovProjection.FunctionProjection needs no modifications to FunctionProjection base class.
local FunctionProjection = FunctionProjectionParent

----------------------------------------------------------------------
-- Vlasov-specific VlasovProjection.MaxwellianProjection extends MaxwellianProjection base class.
local MaxwellianProjection = Proto(MaxwellianProjectionParent)

function MaxwellianProjection:advance(t, inFlds, outFlds)
   local distf = outFlds[1]
   self.project:advance(t, {}, {distf})
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end
end

----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
