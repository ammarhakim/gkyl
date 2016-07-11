-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- One needs to tune loopunroll parameter in the JIT compiler to get
-- some code to run faster. The default is 15. Increasing this can
-- have *dramatic* difference on the run-time. Once an optimum is
-- reached, increasing it further does not help, but may infact make
-- things a little worse. (I am not sure if this should be here. I am
-- assuming that everyone will need to use some updater or the
-- other. Perhaps this should be changed in the main gkyl executable).
jit.opt.start('callunroll=10', 'loopunroll=30')

-- Gkyl modules
local WavePropagation = require "Updater.WavePropagation"
local Bc = require "Updater.Bc"

-- system modules
local xsys = require "xsys"

return xsys.table.union(WavePropagation, Bc)
