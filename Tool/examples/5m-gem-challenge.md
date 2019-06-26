The GEM Reconnection Challenge (GEM Challenge) problem simulated with
the five-moment plasma model. This problem shows that the five-moment
model contains the physics of fast magnetic reconnection. To get the
current sheet structure correct one should include the pressure
tensor. In the five-moment moment model, as the grid is refined, the
current sheet should go to a singular X-point structure and will not
show the characteristic elongation seen in larger system size
simulations.

To plot the out-of-plane electron momentum (or current) do:

```
pgkyl -f 5m-gem_elc_5.bp sel -c3 pl --fix-aspect
```

Note the option `-c3` selects the electron $\rho u_z$. Hence to
compute the current from this one needs to multiply by $q/m$ for the
electrons (in this case $-25$). This can be achieved using:

```
pgkyl -f 5m-gem_elc_5.bp sel -c3 ev "f0 -25 *" pl --fix-aspect
```

This will show the structure of the current sheet at the end of the
simulation. The out-of-plane component of the magnetic field ($B_z$)
can be plotted using:

```
pgkyl -f 5m-gem_field_5.bp  sel -c5 pl --fix-aspect
```

This shows the quadrupolar structure of the Hall magnetic field. Note
that in a ideal MHD simulation this field would be exactly
zero. (Ideal-MHD does not include fast reconnection and one needs at
the minumum Hall currents to get the reconnection rate correct).