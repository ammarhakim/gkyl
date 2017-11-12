# High-precision timers

This module provides high-precision timers (higher than those provided
by the `os` module) using system-specific time APIs.

Load the timer module as:

~~~~~~~ {.lua}
Time = require "Lib.Time"
~~~~~~~

This module provides the following functions:

`Time.time()`
: The wall clock time since some system-specific epoch with about 100
  microsecond precision.

`Time.clock()`
: High precision timer even more accurate than `Time.time()`. It has a
  precision of about 1 microsecond.