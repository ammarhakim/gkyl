-- Gkyl ------------------------------------------------------------------------
--
-- Iterative Poisson solver
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local DataStruct = require "DataStruct"

-- Iterative Poisson solver updater object
local IterPoisson = Proto(UpdaterBase)

-- constructor
function IterPoisson:init(tbl)
   IterPoisson.super.init(self, tbl) -- setup base object

end

return IterPoisson
