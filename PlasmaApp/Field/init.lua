-- for accessing any field object
local FieldBase = require "PlasmaApp.Field.FieldBase"
local MaxwellField = require "PlasmaApp.Field.MaxwellField"
local GkField = require "PlasmaApp.Field.GkField"

return {
  FieldBase = FieldBase.FieldBase,
  FuncFieldBase = FieldBase.FuncFieldBase,
  NoField = FieldBase.NoField,
  MaxwellField = MaxwellField.MaxwellField,
  FuncMaxwellField = MaxwellField.FuncMaxwellField,
  GkField = GkField.GkField,
  --FuncGkField = GkField.FuncGkField,
}
