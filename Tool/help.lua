-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"

-- Usage help for top-level Gkyl executable
local usageString = [[ 

To run a Gkeyll simulation simply specify the name of script you want
to run as the first argument:

  gkyl input.lua

Optionally, you may want to only initialize the simulation without
running it (to plot initial conditions, for example). In this case you
can do:

  gkyl input.lua init

This will simply initialize the simulation but won't run it. To
restart a previously run simulation you should do:

  gkyl input.lua restart

This will restart the simulation from the last stored restart
files. Note that by default Gkeyll writes restart files every 20% of
the total simulation. You can change this using the
'restartFrameEvery' parameter in the top-level App object. By default
this is set to 0.2. To write every 10% you would set this to:

  restartFrameEvery = 0.1,

for example.

Gkeyll provides a set of "tools" to do various things, for example,
show help, run tests or fetch example input files. To run a particular
tool called 'help', for example, specify the name of the tool as the
first argument:

  gkyl help

Some tools may take additional commands or flags. You can find these
by querying their own help system:

  gkyl help -h

 ]]

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Run Gkeyll help system]]
parser:flag("-u --usage", "Print usage information", true)

local args = parser:parse(GKYL_COMMANDS)

-- do stuff
if args.usage then
   io.write (usageString)
end

