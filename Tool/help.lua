-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local KeywordMap = "require Tool.HelpKeywordMap"

-- table of help topics to doc file names
local topics = {
   ['usage'] = 'helpdocs/usage.md',
}

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Run Gkeyll help system]]

parser:argument("optional", "Keyword to search")
   :args("?")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

-- Attempts to open a given URL in the system default browser
local open_cmd
function open_url(url)
   if not open_cmd then
      if package.config:sub(1,1) == '\\' then -- windows
	 open_cmd = function(url)
	    -- Should work on anything since (and including) win'95
	    os.execute(string.format('start "%s"', url))
	 end
      elseif io.popen("uname -s"):read'*a' == "Darwin\n" then -- OSX/Darwin
	 open_cmd = function(url)
	    os.execute(string.format('open "%s"', url))
	 end
      else
	 open_cmd = function(url)
	    os.execute(string.format('xdg-open "%s"', url))
	 end
      end
   end
   open_cmd(url)
end

if not args.optional then
   open_url("https://gkeyll.readthedocs.io/en/latest/")
else
   local kw = args.optional
   print("Looking for keyword: ", kw)
end
