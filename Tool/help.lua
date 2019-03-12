-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local lume = require "Lib.lume"

-- table of help topics to doc file names
local topics = {
   ['usage'] = 'helpdocs/usage.md',
}

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Run Gkeyll help system]]
parser:flag("-v --version", "Gkeyll version information")
parser:flag("-u --usage", "Print usage information")
parser:flag("-l --list", "Show list of help topics")
parser:option("-t --topic", "Show help for topic")
parser:flag("--tools", "Show help for topic")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

-- do stuff
if args.list then
   io.write("Help topics: \n")
   for tn, topic in pairs(topics) do
      io.write(string.format(" %s ", tn))
   end
   io.write("\n")
elseif topics[args.topic] then
   local f = io.open(GKYL_EXEC_PATH .. "/Tool/" .. topics[args.topic], "r")
   local helpStr = f:read("*a")
   io.write (helpStr .. "\n")
elseif args.usage then
   local f = io.open(GKYL_EXEC_PATH .. "/Tool/" .. topics['usage'], "r")
   local helpStr = f:read("*a")
   io.write (helpStr .. "\n")
elseif args.version then
   io.write("Changeset: " .. GKYL_HG_CHANGESET .. "\nBuild date: " .. GKYL_BUILD_DATE .. "\n")
elseif args.tools then
   io.write("Supported tools: \n")
   for tn, tool in pairs(GKYL_TOOLS) do
      io.write(string.format(" %s: %s\n", tn, tool[2]))
   end
else
   io.write(parser:get_help() .. "\n")
end
