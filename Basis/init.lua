-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the basis function module
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local CartModalMaxOrder = require "Basis.CartModalMaxOrder"
local CartModalSerendipity = require "Basis.CartModalSerendipity"

-- system modules
local xsys = require "xsys"

return xsys.table.union(CartModalMaxOrder, CartModalSerendipity)
