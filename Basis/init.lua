-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the basis function module
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local CartModalMaxOrder = require "Basis.CartModalMaxOrder"

-- system modules
local xsys = require "xsys"

return xsys.table.union(CartModalMaxOrder)
