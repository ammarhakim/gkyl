-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the Grid modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CartGrid = require "Grid.CartGrid"

return {
   CartGrid = CartGrid.CartGrid,
   NonUniformCartGrid = CartGrid.NonUniformCartGrid,
}
  


