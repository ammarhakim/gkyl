-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local WavePropagation = require "Updater.WavePropagation"
local Bc = require "Updater.Bc"
local FiveMomentSrc = require "Updater.FiveMomentSrc"

-- system modules
local xsys = require "xsys"

return xsys.table.union(WavePropagation, Bc, FiveMomentSrc)
