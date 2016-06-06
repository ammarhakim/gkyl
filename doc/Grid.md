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
constructor fields. The following methods are provided:

~~~~~~~ {.lua}
grid:ndim()
~~~~~~~

The dimensions of the grid.

~~~~~~~ {.lua}
grid:lower(dir)
grid:upper(dir)
~~~~~~~

The lower/upper bounds in direction `dir`.

~~~~~~~ {.lua}
grid:numCells(dir)
~~~~~~~

The number of cells in direction `dir`.

~~~~~~~ {.lua}
grid:dx(dir)
~~~~~~~

Cell spacing in direction `dir`