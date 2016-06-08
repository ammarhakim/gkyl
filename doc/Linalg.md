# Vectors, matrices and linear algebra

Linear algebra objects and methods can be accessed by loading the
"Linalg" module as follows:

~~~~~~~ {.lua}
Lin = require "Linalg"
~~~~~~~

The following objects are provided b this module.

## `vec`: 1D vector

The `vec` constructor can be used to create vectors of a specifed (but
fixed) length. Vectors are __zero__ indexed, and indexing is performed
using square brackets. For example:

~~~~~~~ {.lua}
v = Lin.vec(3)
for i=0,#v-1 do
  v[i] = (i+0.5)*0.1
end
~~~~~~~

Vectors of fixed-sized C structs can also be created. To do this,
first create a constructor as follows:

~~~~~~~ {.lua}
eulerVec = Lin.new_vec_ct(ffi.typeof("struct {double rho, rhou, E;}"))
~~~~~~~

Now, `eulerVec` is a constructor for creating vectors which store the
struct with `rho`, `rhou` and `E`. For example,

~~~~~~~ {.lua}
v = Lin.vec(3)
for i=0,#v-1 do
  v[i] = (i+0.5)*0.1
end
~~~~~~~

The following methods are provided:

`v:elemType()`
: Returns element LuaJIT FFI type stored in the vector.

`v:data()`
: Returns pxointer to the underlying C data pointer storing the data.

`v:copy()`
: Creates and returns a new vector object with the same contents as `v`.

