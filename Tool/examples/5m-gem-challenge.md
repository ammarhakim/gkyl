The GEM Reconnection Challenge (GEM Challenge) problem simulated with
the five-moment plasma model. This problem shows that the five-moment
model contains the physics of fast magnetic reconnection. To get the
current sheet structure correct one should include the pressure
tensor. See the 10-gem problem for the corresponding input file.

To plot the out-of-plane electron momentum (or current) do:

```
pgkyl -f 5m-gem_elc_5.bp sel -c3 pl --fix-aspect
```

This will show the structure of the current sheet at the end of the
simulation. The out-of-plane component of the magnetic field can be
plotted using:

```
pgkyl -f 5m-gem_field_5.bp  sel -c5 pl --fix-aspect
```

This shows the quadrupolar structure of the Hall magnetic field. Note
that in a ideal MHD simulation this field would be exactly
zero. (Ideal-MHD does not include fast reconnection and one needs at
the minumum Hall currents to get the reconnection rate correct).