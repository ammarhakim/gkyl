# Vectors, matrices and linear algebra

Linear algebra objects and methods can be accessed by loading the
"Linalg" module as follows:

~~~~~~~ {.lua}
Lin = require "Linalg"
~~~~~~~

The following objects are provided b this module.

## `Vec`: 1D vector

The `vec` constructor can be used to create vectors of doubles of a
specifed (but fixed) length. __Vectors are indexed starting at 1__ (as
is standard in Lua) and indexing is performed using square
brackets. For example:

~~~~~~~ {.lua}
v = Lin.Vec(3)
for i=1,#v do
  v[i] = (i-0.5)*0.1
end
~~~~~~~

The pre-defined constructor `IntVec` creates a vector of 32-bit
integers, and `FloatVec` creates a vector of floating point numbers.

Vectors of arbitrary types, including fixed-sized C structs, can also
be created. To do this, first create a constructor as follows:

~~~~~~~ {.lua}
EulerVec = Lin.new_vec_ct(ffi.typeof("struct {double rho, rhou, E;}"))
~~~~~~~

Now, `EulerVec` is a constructor for creating vectors which store the
struct with `rho`, `rhou` and `E`. For example,

~~~~~~~ {.lua}
v = EulerVec(3)
for i=1,#v do
  v[i].rho = 1.0
  v[i].rhou = 2.0
  v[i].E = 3.0
end
~~~~~~~

The following methods are provided:

`v:elemType()`
: Returns element LuaJIT FFI type stored in the vector.

`v:data()`
: Returns pointer to the underlying C data pointer storing the
  data. This should be used with caution. Note that as is standard in
  C the returned data pointer is __indexed starting at 0__.

`v:copy()`
: Creates and returns a new vector object with the same contents as `v`.

