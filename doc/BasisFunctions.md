# Basis functions

Various basis functions are defined in Gkyl can be accessed by loading
the "Basis" modules as follows:

~~~~~~~ {.lua}
Basis = require "Basis" 
~~~~~~~

The basis are defined on a reference cell. For rectangular basis, the
reference cell is a (hyper) cube with side 2, and each side spans
$[-1,1]$.

The following objects are provided by this module.

## `CartModalMaxOrder`: Modal, maximal order basis functions

These basis functions include monomials with a given total order. The
basis are orthonormal. The constructor takes the following parameters:

`ndim`
: Number of dimensions of reference hypercube. Basis up to dimension 6
  are supported.

`polyOrder`
: Polynomial order. This is the total polynomial order of the
  monomials included in the basis set.

The following methods are provided.

`basis:ndim()`
: Number of dimensions of reference element.

`basis:polyOrder()`
: Polynomial order.

`basis:numBasis()`
: Total number of basis functions in the element.

`basis:evalBasis(z, bvalues)`
: The value of the basis functions at reference coordinate `z`. The
  values are returned in the 1-indexed and pre-allocated vector
  `bvalues`.

## `CartModalSerendipity`: Modal, serendipity basis functions

These basis functions include monomials with a specified super-linear
order. The super-linear order of a polynomials is the sum of the
quadratic and greater powers. See. Found. Comput. Math 2011 11:337
Arnold & Awanou. The basis are orthonormal. The constructor takes the
following parameters:

`ndim`
: Number of dimensions of reference hypercube. Basis up to dimension 6
  are supported.

`polyOrder`
: Polynomial order. This is the total polynomial order of the
  monomials included in the basis set.

The following methods are provided.

`basis:ndim()`
: Number of dimensions of reference element.

`basis:polyOrder()`
: Polynomial order.

`basis:numBasis()`
: Total number of basis functions in the element.

`basis:evalBasis(z, bvalues)`
: The value of the basis functions at reference coordinate `z`. The
  values are returned in the 1-indexed and pre-allocated vector
  `bvalues`.

