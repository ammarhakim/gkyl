Two-stream instability of two counter-streaming beams of
electrons. Ions are not evolved and provide a neutralizing background.

You can study the growth of the instability by plotting the $E_x$
component of the electric field energy. To do this:

```
pgkyl -f pgkyl -f two-stream_fieldEnergy_ sel -c0 pl --logy
```

The distribution around the time when the instability goes non-linear
can be plotted using

```
pgkyl -f two-stream_elc_4.bp interp -p 2 -b ms pl
```