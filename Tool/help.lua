-- Gkyl ------------------------------------------------------------------------
--
-- Tool to handle gkyl help system.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"

-- table of help topics to doc file names
local topics = {
   ['usage'] = 'helpdocs/usage.md',
}

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Run Gkeyll help system]]

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

-- Attempts to open a given URL in the system default browser, regardless of Operating System.
local open_cmd -- this needs to stay outside the function, or it'll re-sniff every time...
function open_url(url)
   if not open_cmd then
      if package.config:sub(1,1) == '\\' then -- windows
	 open_cmd = function(url)
	    -- Should work on anything since (and including) win'95
	    os.execute(string.format('start "%s"', url))
	 end
	 -- the only systems left should understand uname...
      elseif io.popen("uname -s"):read'*a' == "Darwin\n" then -- OSX/Darwin
	 open_cmd = function(url)
	    -- I cannot test, but this should work on modern Macs.
	    os.execute(string.format('open "%s"', url))
	 end
      else -- that ought to only leave Linux
	 open_cmd = function(url)
	    -- should work on X-based distros.
	    os.execute(string.format('xdg-open "%s"', url))
	 end
      end
   end
   open_cmd(url)
end

open_url("https://gkeyll.readthedocs.io/en/latest/")
