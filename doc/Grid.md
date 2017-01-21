# Grid objects

Grid objects can be accessed by loading the "Grid" module as follows:

~~~~~~~ {.lua}
Grid = require "Grid"
~~~~~~~  

The following objects are provided by this module.

## `RectCart`: Uniform, cartesian grid

A uniform, cartesian grid can be constructed using the `RectCart`
constructor object. Grid from one to six dimensions can be made. We
refer to the grid dimension as `NDIM`.

For example, the following example creates a 2D grid in parallel on 4
processors:

~~~~~~~ {.lua}
decomp = DecompRegionCalc.CartProd { cuts = {2, 2} }
grid = Grid.RectCart {
   lower = {0.0, 0.0},
   upper = {2.0, 3.0},
   cells = {10, 20},
   decomposition = decomp,
}
~~~~~~~

The grid dimension is determined from the number of entries in the
constructor fields (i.e, `NDIM = #cells`). The grid constructor takes
the following parameters:

`lower`
: Table with coordinates of the lower-left corner of grid.

`upper`
: Table with coordinates of the upper-right corner of grid.

`cells`
: Number of cells in each direction

`periodicDirs` (Optional. Defaults to empty table)
: Periodic directions (indexed from 1 as X, 2 as Y etc).

`decomposition`
: The decomposition object used to decompose grid in parallel. This
  must be specified for use in parallel.

The following methods are provided. Note: __we indicate X-direction
with 1, Y-direction with 2, etc.__

`grid:ndim()`
: The dimensions of the grid.

`grid:comm()`
: Communicator for use in parallel messaging.

`grid:subGridId()`
: For parallel grid, returns integer ID of grid. This number is unique
  for each indpendent portion of the grid, which usually live on
  different procssors.

`grid:decomposedRange()`
: Returns an object of type `DecomposedRange`. This encodes the
  information about the domain decomposition, in particular the
  neighbor information needed in parallel communications.

`grid:lower(dir)`, `grid:upper(dir)`
: The lower/upper bounds in direction `dir`.

`grid:numCells(dir)`
: The number of cells in direction `dir`.

`grid:localRange()`, `grid:globalRange()`
: Local/global index range spanned by grid.

`grid:isDirPeriodic(dir)`
: Returns true if direction `dir` is periodic, false otherwise.

`grid:dx(dir)`
: Cell spacing in direction `dir`.

`grid:cellVolume()`
: Volume of a grid cell. For uniform grid, this is the same for all cells.

## `NonUniformRectCart`: Non-uniform, cartesian grid

A non-uniform cartesian grid can be constructed using the
`NonUniformRectCart` constructor object. A non-uniform grid is
described by the 1D arrays which specify the nodal coordinates in each
direction.

A non-uniform grid can be created in two ways: first, by specifying
functions which map computational space coordinates to physical space
coordinates and second, by directly setting the nodal coordinates
manually. One would use the first method if an analytical mapping
function is known, and the second, if the nodal coordinates are
pre-computed by other means (like reading from a file).

The grid constructor takes the following parameters:

`lower` (Optional. Defaults to  `{0.0, ...}`)
: Table with coordinates of the lower-left corner of grid.

`upper` (Optional. Defaults to `{1.0, ...}`)
: Table with coordinates of the upper-right corner of grid.

`cells`
: Number of cells in each direction

`mappings` (Optional. If missing, nodal spacing is uniform)
: This is a table with `NDIM` number of functions, each of which takes
  a single parameter, `zeta`. The `d`-th function in the list must
  return the physical space coordinate corresponding to the
  computational space coordinate `zeta`.

For example, consider a 1D grid $x\in [0,1]$ with the mapping $x =
\zeta^2$. This gives a grid which is packed near $x=0$, with cell size
getting larger towards $x=1$. To make such a grid one would do:

~~~~~~~ {.lua}
grid = Grid.NonUniformRectCart {
   cells = {16},
   mappings = {
      function (zeta)
        return zeta*zeta
      end,
   }
}
~~~~~~~  

Note that as the `lower` and `upper` fields are missing, the
computational domain is assumed to be $[0,1]$.

The following methods are provided. (Note: __we indicate X-direction
with 1, Y-direction with 2, etc__).

`grid:ndim()`
: The dimensions of the grid.

`grid:lower(dir)`, `grid:upper(dir)`
: The lower/upper bounds in direction `dir`.

`grid:numCells(dir)`
: The number of cells in direction `dir`.

`grid:nodeCoords(dir)`
: Return a 1D vector object (see documentation of 1D vectors for API)
  with the nodal coordinates in direction `dir`. This can be modified
  to manually set the nodal coordinates if needed.

`grid:setIndex(idx)`
: Set index into grid to point to cell `idx`. This method must be
  called before querying the grid for cell size or volume. The index
  can be any object which supports the [] indexing operator. However,
  __please remember that indices start at 1__. This is specially
  important to keep in mind if using a raw ffi object as an index.

`grid:dx(dir)`
: Cell spacing in direction `dir`.

`grid:cellVolume()`
: Volume of a grid cell.


