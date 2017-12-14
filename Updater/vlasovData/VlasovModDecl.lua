-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Vlasov C++ kernel functions based on basis functions
-- and polyOrder
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local VlasovCdef = require "_VlasovCdef"

local _M = {}

function _M.selectVolStream(CDIM, VDIM, basisNm, polyOrder)
end

function _M.selectSurfStream(CDIM, VDIM, basisNm, polyOrder)
end

return _M
