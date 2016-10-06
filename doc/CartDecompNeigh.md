# Computing neighbor relations for parallel communication

To compute the list of neighboring ranks for use in parallel
communication, use the `CartDecompNeigh` object. This object is used
internally and most end-users will not need it. To access this object,
load it as:

~~~~~~~ {.lua}
CartDecompNeigh = require "Lib.CartDecompNeigh"
~~~~~~~

The object is constructed by passing it the decomposition object. See
documentation for `CartDecomp`. Hence, a neighbor calculator object is
created as:

~~~~~~~ {.lua}
local decompNeigh = CartDecompNeigh(decomp)
~~~~~~~

where `decomp` is an already decomposed region.

The following methods are provided.

`decompNeigh:calcFaceCommNeigh(lowerGhost, upperGhost)`
: Compute neighbor information for use in communicating ghost/skin
  cell data with `lowerGhost` number of cells on lower edges and
  `upperGhost` number of cells on upper edges. Note that the set of
  neighbors computed by this method will not include corner neighbors
  and hence should not be used for algorithms which include corner
  transport terms.

`decompNeigh:calcAllCommNeigh(lowerGhost, upperGhost)`
: Compute neighbor information for use in communicating ghost/skin
  cell data with `lowerGhost` number of cells on lower edges and
  `upperGhost` number of cells on upper edges. All neighbors,
  including corner neighbors are included.

`decompNeigh:neighborData(k)`
: Returns table of neighbor sub-domain numbers for sub-domain `k`.


