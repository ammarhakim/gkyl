-- Gkyl ------------------------------------------------------------------------
--
-- Tool to fetch gkyl example input files.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"

-- table of help topics to doc file names
local examples = {
   ['two-stream'] = {
      [[Electrostatic two-stream instability in 1X1V.]],
      'examples/es-two-stream.lua',
      'examples/es-two-stream.md'
   },
   ['5m-gem'] = {
      [[GEM reconnection challenge problem, five-moment model.]],
      'examples/5m-gem-challenge.lua',
      'examples/5m-gem-challenge.md'
   },
}

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("examples")
   :description [[List and fetch Gkeyll examples]]
parser:flag("-l --list", "Show list of help topics")
parser:option("-e --example", "Example to get")
parser:option("-d --describe", "Describe example")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

local function showExamplesList()
   io.write("Available examples: \n")
   for ex, descr in pairs(examples) do
      io.write(string.format(" %s: %s\n", ex, descr[1]))
   end
end

-- do stuff
if args.list then
   showExamplesList()
elseif examples[args.example] then
   local f = io.open(GKYL_EXEC_PATH .. "/Tool/" .. examples[args.example][2], "r")
   local data = f:read("*a")
   --io.write (data .. "\n")
   local f = io.open(args.example .. ".lua", "w")
   f:write (data .. "\n")
elseif examples[args.describe] then
   local f = io.open(GKYL_EXEC_PATH .. "/Tool/" .. examples[args.describe][3], "r")
   local data = f:read("*a")
   io.write (data .. "\n")
else
   showExamplesList()
end
