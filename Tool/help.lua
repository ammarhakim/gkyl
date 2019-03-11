-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"

-- Help doc for top-level Gkyl executable
local usageString = [[
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
files. Note that by default Gkeyll writes restart files every 20% of
the total simulation. You can change this using the
'restartFrameEvery' parameter in the top-level App object. By default
this is set to 0.2. To write every 10% you would set this to:

```lua
  restartFrameEvery = 0.1,
```

for example.

# Using Gkeyll tools

Gkeyll provides a set of "tools" to do various things, for example,
show help, run tests or fetch example input files. To run a particular
tool called 'help', for example, specify the name of the tool as the
first argument:

```
  gkyl help
```

Some tools may take additional commands or flags. You can find these
by querying the tool's own help system:

```
  gkyl help -h
```
 ]]

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Run Gkeyll help system]]
parser:flag("-u --usage", "Print usage information")
parser:flag("-t --tools", "Show list of tools")

local args = parser:parse(GKYL_COMMANDS)

-- function to format text

-- do stuff
if args.usage then
   io.write (usageString)
elseif args.tools then
   io.write("Supported tools: \n")
   for tn, tool in pairs(GKYL_TOOLS) do
      io.write(string.format(" %s: %s\n", tn, tool[2]))
   end
else
   io.write(parser:get_help() .. "\n")
end

