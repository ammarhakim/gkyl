# Data-structures for fields and diagnostics 

Various data-structures defined in Gkyl can be accessed by loading the
"DataStruct" module as follows:

~~~~~~~ {.lua}
DataStruct = require "DataStruct" 
~~~~~~~  

The following objects are provided by this module.

## `Field`: Multi-component fields on Cartesian grids

Gkyl provides a general field object that allows storing multiple
components in each cell of a Cartesian grid. Default data-type that
can be stored in the field is `double`. However, fields of arbitrary
fixed-size types, including struct with named fields, can be
constructed.





