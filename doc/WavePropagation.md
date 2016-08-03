# The finite-volume wave-propagation updater

Finite volume wave propagation updater on rectangular, Cartesian
grid. Non-uniform grids are supported and an embedded BC can also be
included. Before using this updater one must create an equation
object.

An example usage is:

~~~~~~~ {.lua}
local Updater = require "Updater"

eulerEqn = HyperEquation.Euler { gasGamma = 1.4 }
fluidSlvrDir1 = Updater.WavePropagation {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = 0.9,
   updateDirections = {1}
}
~~~~~~~

The constructor takes the following parameters:

`onGrid`
: Grid on which the updater lives

`equation`
: Equation object

`limiter`
: Limiter to use. This is one of "no-limiter", "min-mod", "superbee",
  "van-leer", "monotonized-centered", "beam-warming" or "zero".

`cfl`
: CFL number to use. Depending on the way the updater is used, the CFL
  varies from $1/D$, where $D$ is the dimension of the domain, to
  1.0. The latter is appropriate for the dimensionally split scheme.

`updateDirections`
: Table of directions to update. For example, for a dimensionally
  split scheme one may want to update the X-direction and Y-direction
  separately. Note that $X$ direction is 1, $Y$ is 2, etc.

`cflm` (Optional. Defaults to `1.1*cfl`)
: Maximum CFL to use. This should be set to slightly larger than `cfl`
  to ensure the code does not thrash around as it adjusts the
  time-step.

The object provides a single "advance" method:

`advance(tCurr, dt, inFld, outFld)`
: Advances the solution from `tCurr` to `tCurr+dt`. The `inFld` table
  should have a single input field, holding the solution at
  `tCurr`. The `outFld` table should have a single output field. This
  stores the updated solution. The function returns two values, the
  first a flag indicating if the advanced succeeded and the second a
  suggested time-step for the next call to this method.
