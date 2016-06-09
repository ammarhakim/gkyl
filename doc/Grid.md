# Grid objects

Grid objects can be accessed by loading the "Grid" module as follows:

~~~~~~~ {.lua}
Grid = require "Grid"
~~~~~~~  

The following objects are provided by this module.

## `CartGrid`: Uniform, cartesian grid

A uniform, cartesian grid can be constructed using the `CartGrid`
constructor object. Grid from one to six dimensions can be made. For
example, the following example creates a 2D grid:

~~~~~~~ {.lua}
grid = Grid.CartGrid {
   lower = {0.0, 0.0},
   upper = {2.0, 3.0},
   cells = {10, 20}
}
~~~~~~~

The grid dimension is determined from the number of entries in the
constructor fields. The grid constructor takes the following
parameters:

`lower`
: Table with coordinates of the lower-left corner of grid.

`upper`
: Table with coordinates of the upper-right corner of grid.

`cells`
: Number of cells in each direction

The following methods are provided:

`grid:ndim()`
: The dimensions of the grid.

`grid:lower(dir)`, `grid:upper(dir)`
: The lower/upper bounds in direction `dir`.

`grid:numCells(dir)`
: The number of cells in direction `dir`.

`grid:dx(dir)`
: Cell spacing in direction `dir`.

`grid:cellVolume()`
: Volume of a grid cell. For uniform grid, this is the same for all cells.
