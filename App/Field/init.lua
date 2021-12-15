-- Gkyl ------------------------------------------------------------------------
--
-- App support code: load all field objects
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- for accessing any field object
local FieldBase = require "App.Field.FieldBase"
local AmbipolarSheathField = require "App.Field.AmbipolarSheathField"
local MaxwellField = require "App.Field.MaxwellField"
local GkField = require "App.Field.GkField"
local GkZeroCurrentField = require "App.Field.GkZeroCurrentField"

return {
  FieldBase = FieldBase.FieldBase,
  AmbipolarSheathField = AmbipolarSheathField,
  ExternalFieldBase = FieldBase.ExternalFieldBase,
  FuncFieldBase = FieldBase.ExternalFieldBase, -- for backwards compat
  NoField = FieldBase.NoField,
  MaxwellField = MaxwellField.MaxwellField,
  ExternalMaxwellField = MaxwellField.ExternalMaxwellField,
  FuncMaxwellField = MaxwellField.ExternalMaxwellField, -- for backwards compat
  GkField = GkField.GkField,
  GkGeometry = GkField.GkGeometry,
  GkZeroCurrentField = GkZeroCurrentField,
}
