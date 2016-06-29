# Some notes on implementing hyperbolic equations

To work with the finite-volume scheme each hyperbolic equation object
must implement the following methods. See "Hakim, A., Loverich, J., &
Shumlak, U. (2006). A high resolution wave propagation scheme for
ideal Two-Fluid plasma equations. Journal of Computational Physics,
219(1), 418–442. http://doi.org/10.1016/j.jcp.2006.03.036" and also
LeVeque's book on Finite Volume methods.

`numEquations()`
: Number of equations in system.

`numWaves()`
: Number of waves in system.

`isPositive(q)`
: Check if the input conserved variable vector is in the invariant
  domain of the physical system. (Eg: density and pressure are positive)

`rp(delta, ql, qr, waves, s)`
: Compute the `waves` and speed `s` given the jump `delta`.

`qFluctuations(ql, qr, waves, s, amdq, apdq))`
: Compute the q-fluctuations, given waves and wave speeds. The
  fluctuations in normal cases are defined by $A\Delta Q^- =
  \sum_{s<0} s W$ and $A\Delta Q^+ = \sum_{s>0} s W$. Usually, `qr`
  and `ql` should be ignored. They come into play only when
  implementing Lax or other specific types of fluxes. Shape of
  `waves` matrix is (mwave X meqn).

`rotateToLocalAligned(dir, qIn, qOut)`,`rotateToLocalAligned(dir, qIn, qOut)`
: Rotate input conserved variable vector to a local/global coordinate
  system in direction `dir`.

To ensure an efficient implementation, manually unroll small inner
loops as much as possible. Templates from the `xsys.template` library
can be used to unfold loops automatically.