-- Gkyl ------------------------------------------------------------------------
--
-- Solves the exact Euler Riemann Problem in 1D for ideal gas law
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[
Solves the exact Riemann problem for Euler equations on a 1D
grid. Given the left/right states, number of cells and time, outputs
density, velocity, pressure and internal energy at that time. The
right state can be a vacuum.

You can get a example input file from the tool by doing:

  gkyl exacteulerrp -e > exacteulerinp.lua

(The above command will save the output to the file
"exacteulerinp.lua"). Modify it for your left/right states etc and run
it:

  gkyl exacteulerrp -i exacteulerinp.lua

This will produce the following files:

 exacteulerinp_density.bp, exacteulerinp_velocity.bp,
 exacteulerinp_pressure.bp, exacteulerinp_internalenergy.bp

These can be vizualized/postprocessed as usual with pgkyl.
]]

-- add flags and options
parser:flag("-e --example", "Fetch example input file", false)
parser:option("-i --input", "Input file")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)
