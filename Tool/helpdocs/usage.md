# Running a simulation

To run a Gkeyll simulation simply specify the name of script you want
to run as the first argument (this assumes gkyl executable is in your
path):

```
  gkyl input.lua
```

Optionally, you may want to only initialize the simulation without
running it (to plot initial conditions, for example). In this case you
can do

```
  gkyl input.lua init
```

This will simply initialize the simulation but won't run it. To
restart a previously run simulation you should do:

```
  gkyl input.lua restart
```

This will restart the simulation from the last stored restart
files. Note that by default Gkeyll writes restart files every $1/5$ of
the total simulation. You can change this using the
'restartFrameEvery' parameter in the top-level App object. By default
this is set to 0.2. To write every $1/10$ you would set this to:

```lua
  restartFrameEvery = 0.1,
```

for example.

# Using Gkeyll tools

Gkeyll provides a set of "tools" to perform various actions, for
example, show help, run tests or fetch example input files. To run a
particular tool called 'help', for example, specify the name of the
tool as the first argument:

```
  gkyl help
```

Some tools may take additional commands or flags. You can find these
by querying the tool's own help system:

```
  gkyl help -h
```