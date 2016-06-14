# Range: N-dimensional Cartesian index set

The `Range` object provides a n-dimensional index set and is loaded as
follows:

~~~~~~~ {.lua}
Range = require "Lib.Range"
~~~~~~~

The `Range` constructor takes two tables as input. The first is the
lower indices of the range, and the second the upper indices. Indices
are __inclusive__.

`Range.Range(lower, upper)`
: Create a new range object with start indices given by the `lower` table
  and end indices given by `upper` table. All indices are inclusive.

For example, to create a 2D range object which spans a $10\times 10$
box with start index at $(1,1)$ do:

~~~~~~~ {.lua}
range = Range.Range( {1, 1}, {10, 10} )
~~~~~~~

The following methods are provided.

`range:ndim()`
: The dimensions of range object

`range:lower(dir)`
: Lower index in direction `dir`.

`range:upper(dir)`
: Upper index in direction `dir`. Inclusive.

`range.shape(dir)`
: Number of elements in direction `dir`.

`range.volume()`
: Total number of indices in range.

The range object provides row and column major iterators to loop over
all the indices in the range object. These should be used in a `for`
loop to step over all the indices.

`range:colMajorIter()`, `range:rowMajorIter()`
: Column/row major iterators to loop over all indices in range. The
  indexer returns a `Vec` object with the indices. The return vector
  __should not be modified__.

Example usage for indexer is:

~~~~~~~ {.lua}
for idx in range:colMajorIter() do
  -- do something with idx
end
~~~~~~~



