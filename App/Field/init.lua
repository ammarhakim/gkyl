-- for accessing any field object
local FieldBase = require "App.Field.FieldBase"
local MaxwellField = require "App.Field.MaxwellField"
local GkField = require "App.Field.GkField"

return {
  FieldBase = FieldBase.FieldBase,
  FuncFieldBase = FieldBase.FuncFieldBase,
  NoField = FieldBase.NoField,
  MaxwellField = MaxwellField.MaxwellField,
  FuncMaxwellField = MaxwellField.FuncMaxwellField,
  GkField = GkField.GkField,
  --FuncGkField = GkField.FuncGkField,
}
