# Domain decomposition for Cartesian grids

Objects that decompose rectangular grids into smaller domains can be
accessed by loading the `CartDecomp` module:

~~~~~~~ {.lua}
DecompRegionCalc = require "Lib.CartDecomp"
~~~~~~~

This module provides the following objects for decomposing Cartesian
grids using different methods.

## `CartProd`: Decomposition specified by "cuts"

The `CartProd` object allows decomposing a Cartesian grid by
specifying the number of blocks to use in each direction. For example,
to split a 2D grid into $4\times 8$ blocks (i.e. total of 32
sub-domains) one would do:

~~~~~~~ {.lua}
decomp = DecompRegionCalc.CartProd { cuts = {4, 8} }
~~~~~~~

Note that the dimension of the region to decompose is determined from
the size of the `cuts` parameter. The constructor takes no other
parameters.

The following methods are provided.

`decomp:ndim()`
: Dimension of decomposition.

`decomp:cuts(dir)`
: Cuts in direction `dir`.

`decomp:numSubDomains()`
: Total number of sub-domains in decomposition.
  
`decomp:decompose(range)`
: Decomposes `range` into smaller regions. Returns `true` if
  decomposition was successful, `false` otherwise.

After the decompose method has been called, the following methods can
be used to access the decomposition.

`decomp:subDomain(k)`
: The k-th sub-domain in the system, indexed from 1.

`decomp:boundarySubDomainIds(dir)`
: Boundary sub-domain IDs along direction `dir`. The returned Ids are
  stored as a (one-indexed) vector (Linalg.Vec) of pairs. The first
  element in the pair is called `lower` and the second element is
  called `upper`. This function is mostly useful to determine
  communication neighbors in periodic directions.

