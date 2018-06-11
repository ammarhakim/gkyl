-- Gkyl ------------------------------------------------------------------------
--
-- Mapped Cartesian grids: computational space is rectangular but
-- physical space need not be.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"

-- Gkyl libraries
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local RectCart = require "Grid.RectCart"

-- MappedCartGrid --------------------------------------------------------------
--
-- Grid determined by a coordinate mapping between uniform rectangular
-- computational space to physical space.
--------------------------------------------------------------------------------

local MappedCart = Proto(RectCart) -- extends RectCart
-- ctor
function MappedCart:init(tbl)
   MappedCart.super.init(self, tbl)
end

return MappedCart
