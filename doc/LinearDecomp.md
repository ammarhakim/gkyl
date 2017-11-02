# Decomposition for a 1D domain

Objects that decompose 1D domains into smaller domains equitably can
be accessed by loading the `LinearDecomp` module:

~~~~~~~ {.lua}
local LinearDecomp = require "Lib.LinearDecomp"
~~~~~~~

This module provides the following objects for decomposing 1D
domains.

## `LinearDecomp`: Decomposition a domain into number of sub-domains

The `LinearDecomp` object allows decomposing a 1D domain
equitably. For example

~~~~~~~ {.lua}
decomp = LinearDecomp.LinearDecomp { domSize = 100, numSplit = 10 }
~~~~~~~

will divide 100 elements into 10 pieces. The object ensures that every
sub-domain is at most off by 1 from each other.

The constructor takes the following parameters:

`domSize`
: Number of elements to split.

`numSplit` (Optional. Defaults to false)
: Number of splits to make

The following methods are provided.

`decomp:domSize()`
: Size of domain.

`decomp:numSplit()`
: Number of splits.

`decomp:lower(n)`
: Lower index of the n-th sub-domain.

`decomp:shape(n)`
: Number of elements in the n-th sub-domain.
  
`decomp:upper(n)`
: Upper index of the n-th sub-domain.

## `LinearDecompRange`: Decomposition a range into number of "linear" sub-domains

The `LinearDecompRange` object allows decomposing a range object into
linear sub-domain. It calculates the start index for each sub-domain
and its size.

~~~~~~~ {.lua}
r = Range.Range({1, 1, 1}, {100, 100, 100})
decomp = LinearDecomp.LinearDecompRang { range = r, numSplit = 10 }
~~~~~~~

will split a $100^3$ range into 10 linear pieces.

The constructor takes the following parameters:

`range`
: The range to split.

`numSplit` (Optional. Defaults to false)
: Number of splits to make

The following methods are provided.

`decomp:domSize()`
: Size of domain.

`decomp:numSplit()`
: Number of splits.

`decomp:shape(n)`
: Number of elements in the n-th sub-domain.
  
`decomp:rowStartIndex(n)`
: Start index for row major indexing for n-th sub-domain.

`decomp:colStartIndex(n)`
: Start index for column major indexing for n-th sub-domain.