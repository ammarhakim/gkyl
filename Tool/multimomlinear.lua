-- Gkyl ------------------------------------------------------------------------
--
-- Linear dispersion solver for multi-moment multifluid equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local DataStruct = require "DataStruct"
local lfs = require "Lib.lfs"
local ffi = require "ffi"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("multimomlinear")
   :description [[

Solves the linear dispersion relation for multi-moment multifluid
equations. Arbitrary number of species are supported, and the field
equations can either be full Maxwell equations or Poisson equations
for electrostatic problems. Each plasma species can be either an
isothermal fluid (including cold fluid), five-moment or a ten-moment
fluid. Some support for Hammett-Perkin type collisionless closures is
provided.

For full details on how the dispersion solver works see Developer
Documentation on Gkeyll RTD website.

]]
