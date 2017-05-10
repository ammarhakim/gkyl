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

`syncCorners` (Options. Default is `false`)
: If set to true, corner ghost cells are also synchornized on a
 `sync()` call. For most algorithms, this flag should be set to
 `false`, unless the stencil includes use of corner cells.

`syncPeriodicDirs` (Options. Default is `true`)
: If set to true, periodic directions (specified in Grid object) are
  also synchornized on a `sync()` call. In almost all cases the
  default (`true`) is the correction option. In some rare cases one
  may want a field that lives on a periodic domain but yet skip the
  periodic boundary conditions for that field.

The following methods are provided.

`field:elemType()`
: Returns element LuaJIT FFI type stored in the field.

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

`field:genIndexer()`
: A linear indexer object that allows converting a N-dimensional index
  into an integer. The indexer returned by this method takes a vector,
  unliked the indexer returned by the `indexer()` method, which takes
  (i,j,..) index. This is used to access elements in the field.

`field:copy(fIn)`
: Copies data from `fIn` into `field`.

`field:get(k)`
: Get a pointer to field components at location `k`. The value of `k`
  must be determined by the indexer returned by the `indexer()`
  object. The returned object can be indexed as a vector. The
  `numComponents` field can be used to determine the number of
  components in the vector.

`field:getDataPtrAt(k)`
: Get a raw pointer to field components at location `k`. The value of
  `k` must be determined by the indexer returned by the `indexer()`
  object. Note this returns a pointer to the underlying raw memory
  location and hence must by indexed starting with 0. In general, this
  is not a safe method to use, but is provide to make interacting with
  C easier.

`field:fill(k, ptr)`
: Sets a pointer to field components at location `k`. The `ptr` must
  be created by a previous `field:get(0)` call, and the value of `k`
  must be determined by the indexer returned by the `indexer()`
  object. This method is useful in inner loops were using `get` method
  lead to memory fragmentation.

`field:sync()`
: Synchornize values in ghost cells by copy data from neighboring
  ranks' skin cells.

`field:write(outNm, tmStamp)`
: Write data in field to ADIO BP file `outNm`. The parameter `tmStamp`
  is simulation time at which data is written.

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

We can also access the field data in a dimensionally independent
manner. In this approach the above example can be written as:

~~~~~~~ {.lua}
local localRange = field:localRegion()
local indexer = field:genIndexer()
for idx in localRange:colMajorIter() do
   local fitr = field:get(indexer(idx))
   fitr[1] = 1
   fitr[2] = 2
   fitr[3] = 3
end
~~~~~~~

Note the use of the `genIndexer()` method to get the dimensionally
independent indexer.

The field can be used to store data of arbitrary types, including
fixed-size C structs. To do this, first create a new field constructor
as follows:

~~~~~~~ {.lua}
EulerField = DataStruct.new_field_ct(ffi.typeof("struct {double rho, rhou, E;}"))
~~~~~~~

Now, using this, a field can be created:

~~~~~~~ {.lua}
field = EulerField {
   onGrid = grid,
   ghost = {1, 1},
}
~~~~~~~

To loop over the field one can do, for example, using the
dimensionally independent technique:

~~~~~~~ {.lua}
local localRange = field:localRegion()
local indexer = field:genIndexer()
for idx in localRange:colMajorIter() do
   local fitr = field:get(indexer(idx))
   fitr[1].rho = 1
   fitr[1].rhou = 0
   fitr[1].E = 3
end
~~~~~~~

## `DynVector`: Dynamically adjustable 1D array

Gkyl provides a dynamic 1D array to store small amounts of data. The
usual application of such `DynVector` objects is to store
time-dependent diagnostic information like energy history, field
values in a cell, etc. As such, a `DynVector` allows storing
diagnostic at much higher frequency than the I/O frequency of the
simulation.

A `DynVector` takes a single parameter, `numComponents` in its
constructor:

~~~~~~~ {.lua}

emEnergy = DataStruct.DynVector { numComponents = 1 }
~~~~~~~

The following methods are provided.

`dynVec:numComponents()`:
: Number of components.

`dynVec:appendData(tm, vals)`:
: Append data recorded at time `tm` to the end of `dynVec`. The data
  to append must be provided in `vals`, which is a 1-indexed array (or
  table).

`dynVec:removeLast()`:
: Remove the last time and value added to `dynVec`. The removed values
  are returned as a time, value pair.

`dynVec:lastTime()`:
: Return the last time a value was inserted.

`dynVec:lastData()`:
: Returns a pair of values, the first the time and the second the
  values last inserted into the `dynVec`.

`dynVec:timeMesh()`
: Returns the complete 1D array with times at which data was inserted
  into the `dynVec`.

`dynVec:data()`
: Returns the complete 1D array with values inserted into the
  `dynVec`. The returned array is 1-indexed and the [] operator
  returns a 1-indexed array of size `numComponents`.

`dynVec:clear()`
: Clears data stored in the `dynVec`.

`dynVec:write(outNm, tmStamp)`
: Write the data to ADIO BP file `outNm`. The file is time-stamped
  with `tmStamp`. Note that a write clears out the data in `dynVec`,
  hence on subsequent writes only the data stored since the last write
  is written out.
