-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the DataStruct modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CartField = require "DataStruct.CartField"
local DynVector = require "DataStruct.DynVector"
local SparseTriples = require "DataStruct.SparseTriples"

return {
   new_field_ct = CartField.new_field_ct,
   Field = CartField.Field,
   DynVector = DynVector.DynVector,
   SparseTriples = SparseTriples.SparseTriples
}


