# Vectors, matrices and linear algebra

Linear algebra objects and methods can be accessed by loading the
"Linalg" module as follows:

~~~~~~~ {.lua}
Lin = require "Linalg"
~~~~~~~

The following objects are provided b this module.

## `vec`: 1D vector

The `vec` constructor can be used to create vectors of a specifed
length:

~~~~~~~ {.lua}
v = vec(n)
~~~~~~~

This will create a vector of `n` doubles.

