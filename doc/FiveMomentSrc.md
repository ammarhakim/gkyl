# Updating five-moment source terms

The `FiveMomentSrc` updater can be used to update five-moment source
terms. The terms can be updated using an explicit scheme or an
implicit scheme. For most problems its best to use an implicit scheme,
and that is the default behavior.

An example usage is:

~~~~~~~ {.lua}
local Updater = require "Updater"

srcUpdater = Updater.FiveMomentSrc {
   onGrid = grid,
   numFluids = 2,
   charge = {-1.0, 1.0},
   mass = {1.0, 1/1836.2},
   epsilon0 = 1.0,
   scheme = "implicit", -- one of "implicit" or "explicit"
   hasStaticField = false, -- do we have static EM field?
   gravity = 0.0, -- gravitational force
   dir = 1, -- direction of force
}
~~~~~~~

The constructor takes the following parameters:

`onGrid`
: Grid on which the updater lives

`numFluids`
: The number of fluids to update

`charge`
: Table of size `numFluid` with species charges

`mass`
: Table of size `numFluid` with species mass

`epsilon0`
: Premittivity of free space

`scheme` (Optional. Defaults to "time-centered")
: Scheme to use. One of "ssp-rk3", "modified-boris", "backward-euler"  or "time-centered"

`hasStaticField` (Optional. Defaults to `false`)
: If a static electric and magnetic field should be included. The
  field itself, if present, should be provided as an input field. See
  below.

`gravity` (Optional. Defaults to 0.0)
: Gravitational acceleration, if any

`dir` (Optional. Defaults to 1, i.e. X-direction.)
: Direction of gravitational acceleration. X-direction is 1,
  Y-direction 2 and Z-direction 3.

`elcErrorSpeedFactor` (Optional. Defaults to 0.0)
: Factor for propagation of error in electric field divergence

`mgnErrorSpeedFactor` (Optional. Defaults to 0.0)
: Factor for propagation of error in magnetic field divergence

`hasPressure` (Optional. Defaults to `true` for each fluid)
: Table of size `numFluid`, with each entry indicating if the
  corresponding fluid is isothermal, and hence has no pressure
  equation.
