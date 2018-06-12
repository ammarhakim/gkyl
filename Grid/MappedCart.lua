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
local diff = require "sci.diff"

-- MappedCartGrid --------------------------------------------------------------
--
-- Grid determined by a coordinate mapping between uniform rectangular
-- computational space to physical space.
--------------------------------------------------------------------------------

local MappedCart = Proto(RectCart) -- extends RectCart
-- ctor
function MappedCart:init(tbl)
   MappedCart.super.init(self, tbl)

   -- store pointer to mapping function
   self._mapc2p = tbl.mapc2p

   local xp = Lin.Vec(3)
   -- analytical differentiation functions for mapc2p
   self._coordDiff = {
      diff.gradientf(
	 function (xc)
	    local x1, x2, x3 = self._mapc2p(xc)
	    return x1
	 end,
	 self:ndim()
      ),
      diff.gradientf(
	 function (xc)
	    local x1, x2, x3 = self._mapc2p(xc)
	    return x2
	 end,
	 self:ndim()
      ),
      diff.gradientf(
	 function (xc)
	    local x1, x2, x3 = self._mapc2p(xc)
	    return x3
	 end,
	 self:ndim()
      ),      
   }

   -- stuff for use in various methods
   self._xc = Lin.Vec(self:ndim())
end

function MappedCart:mapc2p(xc)
   return self._mapc2p(xc)
end

-- internal methods, not to be used directly by user
function MappedCart:_calcMetric_1d(xc)
   local d1 = Lin.Vec(1)
   self._coordDiff[1](xc, d1)

   local g = Lin.Vec(1)
   g[1] = d1[1]^2
   return g
end
function MappedCart:_calcMetric_2d(xc)
   local d1, d2 = Lin.Vec(2), Lin.Vec(2)
   self._coordDiff[1](xc, d1)
   self._coordDiff[2](xc, d2)

   local g = Lin.Vec(3)
   g[1] = d1[1]^2 + d2[1]^2 -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] -- g_12 = g_21
   g[3] = d1[2]^2 + d2[2]^2 -- g_22

   return g
end
function MappedCart:_calcMetric_3d(xc)
   local d1, d2, d3 = Lin.Vec(3), Lin.Vec(3), Lin.Vec(3)
   self._coordDiff[1](xc, d1)
   self._coordDiff[2](xc, d2)
   self._coordDiff[3](xc, d3)

   local g = Lin.Vec(6)
   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2 -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] + d3[1]*d3[2]  -- g_12 = g_21
   g[3] = d1[1]*d1[3] + d2[1]*d2[3] + d3[1]*d3[3]  -- g_13 = g_31
   g[4] = d1[2]^2 + d2[2]^2 + d3[2]^2 -- g_22
   g[5] = d1[2]*d1[3] + d2[2]*d2[3] + d3[2]*d3[3]  -- g_23 = g_32
   g[6] = d1[3]^2 + d2[2]^2 + d3[3]^2 -- g_33

   return g   
end

-- Computes metric tensor
function MappedCart:calcMetric(xc)
   local ndim = self:ndim()
   if ndim == 1 then
      return self:_calcMetric_1d(xc)
   elseif ndim == 2 then
      return self:_calcMetric_2d(xc)
   elseif ndim == 3 then
      return self:_calcMetric_3d(xc)
   else
      assert(false, "MappedCart does not support more than 3 dimensions!")
   end
end

function MappedCart:cellCenter(xp)
   self.super.cellCenter(self, self._xc)
   if self:ndim() == 1 then
      xp[1] = self._mapc2p(self._xc)
   elseif self:ndim() == 2 then
      xp[1], xp[2] = self._mapc2p(self._xc)
   elseif self:ndim() == 3 then
      xp[1], xp[2], xp[3] = self._mapc2p(self._xc)
   end
end

function MappedCart:write(fName)
   
end

return MappedCart
