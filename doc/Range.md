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

`range:lowerAsVec()`
: Lower indices returned as a NDIM sized vector.

`range:upperAsVec()`
: Upper indices returned as a NDIM sized vector. Inclusive.

`range.shape(dir)`
: Number of elements in direction `dir`.

`range.volume()`
: Total number of indices in range.

`range.copy(rng)`
: Copy range specified by `rng` into this range, i.e `range`.

`range:shorten(dir)`
: Reduces the upper limit in direction `dir` such that the shape in
  that direction is 1. This allows looping over all indices orthogonal
  to `dir`.

`range:extend(loExt, upExt)`
: Returns a new Range object with additional `loExt` indices added
  to all lower edges and `upExt` added to all upper edges.

`range:extend(dir, loExt, upExt)`
: Returns a new Range object with additional `loExt` indices added to
  lower edge and upExt` added to upper edge in direction, `dir`.

`range:shift(offsets)`
: Returns a new Range object shifted by `offsets`

`range:shiftInDir(dir, offset)`
: Returns a new Range object shifted by `offset` in direction `dir`.

`range:lowerSkin(dir, nGhost)/upperSkin(dir, nGhost)`
: Returns a new Range object representing the lower/upper skin
  regions in direction `dir`. The width of the skin region is
  `nGhost`.

`range:lowerGhost(dir, nGhost)/upperGhost(dir, nGhost)`
: Returns a new Range object representing the lower/upper ghost
  regions in direction `dir`. The width of the ghost region is
  `nGhost`.

`range:intersect(rgn)`
: Returns intersection of `range` with `rgn`. If the intersection is
  empty, the lower/upper are such that the volume of the region is
  zero.

`range:isIntersectionEmpty(rgn)`
: Returns true if intersection of `range` and `rgn` is empty, and
  false otherwise.

`range1 == range2`
: Return `true` if ranges are the same, false otherwise.

The range object provides row and column major iterators to loop over
all the indices in the range object. These should be used in a `for`
loop to step over all the indices.

`range:colMajorIter()`, `range:rowMajorIter()`
: Column/row major iterators to loop over all indices in range. The
  indexer returns a `Vec` object with the indices. The return vector
  __should not be modified__.

Example usage for iterators is:

~~~~~~~ {.lua}
for idx in range:colMajorIter() do
  -- do something with idx
end
~~~~~~~

The `Range` module also provides two sets of methods to map a
N-dimensional index to a single number than can be used to index a
flat array. These methods take a range and return a function which
maps indices to a number:

`Range.makeRowMajorIndexer(range)`, `Range.makeColMajorIndexer(range)`
: Returns a function that takes a N-dimensional $(i,j,...)$ index and
  returns a integer. The returned integer lies in the range $[1, V]$,
  where $V$ is the volume of the supplied range.

`Range.makeRowMajorGenIndexer(range)`, `Range.makeColMajorGenIndexer(range)`
: Returns a function that takes a N-dimensional index, $I$ and returns
  a integer. The returned integer lies in the range $[1, V]$, where
  $V$ is the volume of the supplied range.

Example usage of the first set of functions

~~~~~~~ {.lua}
local range = Range.Range({1, 1}, {10, 5})
local indexer = Range.makeRowMajorIndexer(range)
for i = range:lower(1), range:upper(1) do
   for j = range:lower(2), range:upper(2) do
      local linIndex = indexer(i, j)
      -- do something with linIndex
   end
end
~~~~~~~

Example usage of the second set of functions

~~~~~~~ {.lua}
local range = Range.Range({1, 1}, {10, 5})
local indexer = Range.makeRowMajorGenIndexer(range)
for idx in range:rowMajorIter() do
   local linIndex = indexer(idx)
   -- do something with linIndex
end
~~~~~~~