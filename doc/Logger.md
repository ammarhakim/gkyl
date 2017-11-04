# Simulation message logging

This module provides an object to log messages from a
simulation. Messages are written to console and, optionally, to a log
file. When running in parallel, a logging rank can be specified.

Load the logger module as:

~~~~~~~ {.lua}
local Logger = require "Lib.Logger"
~~~~~~~

To create a new logger do, for example,

~~~~~~~ {.lua}
log = Logger { logToFile = true }
~~~~~~~

This will create a logger that will write to a log file in addition to
the console. The logger constructor takes the following parameters:

`writeRank` (Optional. Defaults to 0)
: MPI rank to use for logging

`comm` (Optional. Defaults to Mpi.COMM_WORLD)
: MPI communicator on which rank is being chosen.

`logToFile` (Optional. Defaults to false)
: If set to true, a log file with extention .log will be produced.

The logger takes a single string message to log:

~~~~~~~ {.lua}
log("This is a log message")
~~~~~~~

Use the string.format method (or string concatenation) to produce a
formatted message.

