-- for accessing any field object
local FieldBase = require "GkylApp.Field.FieldBase"
local MaxwellField = require "GkylApp.Field.MaxwellField"
local GkField = require "GkylApp.Field.GkField"

return {
  FieldBase = FieldBase.FieldBase,
  FuncFieldBase = FieldBase.FuncFieldBase,
  NoField = FieldBase.NoField,
  MaxwellField = MaxwellField.MaxwellField,
  FuncMaxwellField = MaxwellField.FuncMaxwellField,
  GkField = GkField.GkField,
  --FuncGkField = GkField.FuncGkField,
}
