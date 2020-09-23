-- Gkyl ------------------------------------------------------------------------
--
-- App support code: load all field objects
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- for accessing any field object
local FieldBase = require "App.Field.FieldBase"
local MaxwellField = require "App.Field.MaxwellField"
local GkField = require "App.Field.GkField"

return {
  FieldBase = FieldBase.FieldBase,
  ExternalFieldBase = FieldBase.ExternalFieldBase,
  FuncFieldBase = FieldBase.ExternalFieldBase, -- for backwards compat
  NoField = FieldBase.NoField,
  MaxwellField = MaxwellField.MaxwellField,
  ExternalMaxwellField = MaxwellField.ExternalMaxwellField,
  FuncMaxwellField = MaxwellField.ExternalMaxwellField, -- for backwards compat
  GkField = GkField.GkField,
  GkGeometry = GkField.GkGeometry,
}
