# (Perfectly Hyperbolic) Maxwell equations for electromagnetism

The `PhMaxwell` object can be used to solve Maxwell equations in their
"perfectly hyperbolic" form. These equations take into account the
divergence relations, and "clean" divergence errors by advecting them
out of the domain. For theory behind these equations and their
properties see, for example: C.-D Munz, P. Omnes, R. Schneider and
E. Sonnendruer and U. Voss, "Divergence Correction Techniques for
Maxwell Solvers Based n a Hyperbolic Model", Journal of Computational
Physics, 161, 484-511, 2000.

The `PhMaxwell` object is in the `Eq` module, and can be loaded as

~~~~~~~ {.lua}
HyperEquation = require "Eq"
~~~~~~~

To construct the object, for example, do:

~~~~~~~ {.lua}
maxwell = HyperEquation.PhMaxwell {
   lightSpeed = 299792458,
   elcErrorSpeedFactor = 1.0,
   mgnErrorSpeedFactor = 1.0
}   
~~~~~~~

The constructor takes the following parameters:

`lightSpeed`
: Speed of light. This need not be in any specific unit, but care must
  be taken to ensure consistency with other parameters in the
  simulation.

`elcErrorSpeedFactor` (Optional. Defaults to 0.0)
: Speed-factor at which divergence errors in electric field are
  advected out of the domain. The advection speed is
  elcErrorSpeedFactor*lightSpeed. A good choice is set this to 1.0.

`mgnErrorSpeedFactor` (Optional. Defaults to 0.0)
: Speed-factor at which divergence errors in magnetic field are
  advected out of the domain. The advection speed is
  `mgnErrorSpeedFactor*lightSpeed`. A good choice is set this to 1.0.

Note: The perfectly hyperbolic Maxwell equations often do not work
very well when `elcErrorSpeedFactor` is set. It is sometimes better
not to correct the electric field divergence error.

The following methods are provided:

`maxwell:numEquations()`
: Number of equations in system. Is always 8.

`maxwell:numWaves()`
: Number of waves in system. Is always 6.

`maxwell:flux(dir, qIn, fOut)`
: Compute the flux (output in `fOut`) given the conserved variable `qIn`.

`maxwell:rp(dir, delta, ql, qr, waves, s)`
: Computes the `waves` and speed `s` given the jump `delta`. Shape of
  `waves` matrix is (mwave X meqn).

`maxwell:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)`
: Computes the q-fluctuations, given waves and wave speeds. The
  fluctuations in normal cases are defined by $A\Delta Q^- =
  \sum_{s<0} s W$ and $A\Delta Q^+ = \sum_{s>0} s W$. Shape of `waves`
  matrix is (mwave X meqn).