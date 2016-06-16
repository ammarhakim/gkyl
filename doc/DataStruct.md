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

A field lives on a grid. Hence, to create a field, first you must
create the grid on which it lives. In general, several fields will
live on the same grid, and most simulations will have a single grid
and several fields. For example, to create a field storing three
components (say the components of the electric field) on a 2D grid,
one would do:

~~~~~~~ {.lua}
grid = Grid.RectCart {
   lower = {0.0, 0.0},
   upper = {1.0, 1.0},
   cells = {32, 32}
}

field = DataStruct.Field {
   onGrid = grid,
   numComponents = 3,
   ghost = {1, 1},
}
~~~~~~~

This will create a $32x32$ field with 3 components with a single layer
of ghost cells on each face of the grid.

The field constructor takes the following parameters:

`onGrid`
: Grid on which the field lives.

`numComponents` (Optional. Defaults to  1)
: Number of components in each grid cell.

`ghost` (Optional. Defaults to  `{0, 0}`)
: Number of ghost cells on the "left" and "right" faces in each direction.

`layout` (Optional. Defaults to "col-major")
: Layout of the data. The default is column-major order. Only other option
  is "row-major". 

The following methods are provided.

`field:ndim()`
: Field dimension. This is the same as the grid dimension.

`field:numComponents()`
: Number of components in each cell in grid.

`field:layout()`
: Layout of field. Either "row-major" or "col-major".

`field:lowerGhost()`, `field:upperGhost()`
: Number of ghost cells on the lower/upper face of grid.

`field:localRange()`, `field:globalRange()`
: Local/global index range spanned by field.

`field:localExtRange()`, `field:globalExtRange()`
: Extended local/global index range spanned by field. This includes
  the indices of the ghost cells.

`field:size()`
: The total number of elements stored in field, including ghost cells.

`field:indexer()`
: A linear indexer object that allows converting a N-dimensional index
  into an integer. This is used to access elements in the field.

`field:get(k)`
: Get the field components at location `k`. The value of `k` must be
  determined by the indexer returned by the `indexer()` object.

To illustrate the use of the `indexer()` and `get()` methods to access
elements in the grid, consider the following code:

~~~~~~~ {.lua}
local localRange = field:localRegion()
local indexer = field:indexer()
for i = localRange:lower(1), localRange:upper(1) do
   for j = localRange:lower(2), localRange:upper(2) do
      local fitr = field:get(indexer(i,j))
      fitr[1] = i+2*j+1
      fitr[2] = i+2*j+2
      fitr[3] = i+2*j+3
   end
end
~~~~~~~

Note the use of the `indexer()` method to `get()` access to the data in
the (i,j) cell. Once the data is fetched, the `fitr` can be indexed to
get the components stored in that cell.