-- Gkyl ------------------------------------------------------------------------
--
-- Mapped Cartesian grids: computational space is rectangular but
-- physical space need not be.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
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

   -- function that maps the computational coordinates to cartesian coordinates
   -- of the form 
   -- X, Y, Z = f({xcomp1, xcomp2, xcomp3})
   self._mapc2p = tbl.mapc2p

   -- determine how many values mapc2p returns
   self._rdim = #{ self._mapc2p(self._lower) }

   -- analytical differentiation functions for mapc2p
   self._coordDiff = {
      diff.gradientf(
	 function (xc)
	    local x1, _, _ = self._mapc2p(xc)
	    return x1
	 end,
	 self:ndim()
      ),
      diff.gradientf(
	 function (xc)
	    local _, x2, _ = self._mapc2p(xc)
	    return x2
	 end,
	 self:ndim()
      ),
      diff.gradientf(
	 function (xc)
	    local _, _, x3 = self._mapc2p(xc)
	    return x3
	 end,
	 self:ndim()
      ),
   }

   -- stuff for use in various methods
   self._xc = Lin.Vec(self:ndim())
   self._d1, self._d2, self._d3 = Lin.Vec(self:ndim()), Lin.Vec(self:ndim()), Lin.Vec(self:ndim())
end

--function MappedCart:id() return "mapped" end
function MappedCart:rdim() return self._rdim end

function MappedCart:mapc2p(xc)
   return self._mapc2p(xc)
end

-- internal methods, not to be used directly by user
function MappedCart:_calcMetric_1d(xc, g)
   local d1 = self._d1
   self._coordDiff[1](xc, d1)

   g[1] = d1[1]^2
end
function MappedCart:_calcMetric_2d_r2(xc, g)
   local d1, d2 = self._d1, self._d2
   self._coordDiff[1](xc, d1)
   self._coordDiff[2](xc, d2)

   g[1] = d1[1]^2 + d2[1]^2 -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] -- g_12 = g_21
   g[3] = d1[2]^2 + d2[2]^2 -- g_22
end
function MappedCart:_calcMetric_2d_r3(xc, g)
   local d1, d2, d3 = self._d1, self._d2, self._d3
   self._coordDiff[1](xc, d1)
   self._coordDiff[2](xc, d2)
   self._coordDiff[3](xc, d3)

   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2 -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] + d3[1]*d3[2] -- g_12 = g_21
   g[3] = d1[2]^2 + d2[2]^2 + d3[2]^2 -- g_22
end
function MappedCart:_calcMetric_3d(xc, g)
   local d1, d2, d3 = self._d1, self._d2, self._d3
   self._coordDiff[1](xc, d1)
   self._coordDiff[2](xc, d2)
   self._coordDiff[3](xc, d3)

   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2 -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] + d3[1]*d3[2]  -- g_12 = g_21
   g[3] = d1[1]*d1[3] + d2[1]*d2[3] + d3[1]*d3[3]  -- g_13 = g_31
   g[4] = d1[2]^2 + d2[2]^2 + d3[2]^2 -- g_22
   g[5] = d1[2]*d1[3] + d2[2]*d2[3] + d3[2]*d3[3]  -- g_23 = g_32
   g[6] = d1[3]^2 + d2[3]^2 + d3[3]^2 -- g_33
end

-- Computes (covariant) metric tensor g_ij
function MappedCart:calcMetric(xc, gOut)
   local ndim = self:ndim()
   if ndim == 1 then
      self:_calcMetric_1d(xc, gOut)
   elseif ndim == 2 then
      if self._rdim == 2 then
	 self:_calcMetric_2d_r2(xc, gOut)
      else
	 self:_calcMetric_2d_r3(xc, gOut)
      end
   elseif ndim == 3 then
      self:_calcMetric_3d(xc, gOut)
   else
      assert(false, "MappedCart does not support more than 3 dimensions!")
   end
end

-- Compute jacobian = det(g_ij)^(1/2)
function MappedCart:calcJacobian(xc)
   local g = {}
   local ndim = self:ndim()
   self:calcMetric(xc, g)
   local jacobian
   if ndim == 1 then 
      jacobian = 1
   elseif ndim == 2 then 
      jacobian = math.sqrt(g[1]*g[3] - g[2]*g[2]) 
   elseif ndim == 3 then
      jacobian = math.sqrt(-g[3]^2*g[4] + 2*g[2]*g[3]*g[5] - g[1]*g[5]^2 - g[2]^2*g[6] + g[1]*g[4]*g[6])
   end
   return jacobian
end  

-- Computes contravariant metric tensor g^ij = (g_ij)^-1
function MappedCart:calcContraMetric(xc, gContraOut)
   local ndim = self:ndim()
   local g = {}
   self:calcMetric(xc, g)
   if ndim == 1 then 
      gContraOut[1] = 1/g[1]
   elseif ndim == 2 then
      local det = self:calcJacobian(xc)^2
      gContraOut[1] = g[3]/det   -- g^11
      gContraOut[2] = -g[2]/det  -- g^12
      gContraOut[3] = g[1]/det   -- g^22
   elseif ndim == 3 then
      local det = self:calcJacobian(xc)^2
      gContraOut[1] = (g[4]*g[6]-g[5]^2)/det     -- g^11
      gContraOut[2] = (g[3]*g[5]-g[2]*g[6])/det  -- g^12 = g^21
      gContraOut[3] = (g[2]*g[5]-g[3]*g[4])/det  -- g^13 = g^31
      gContraOut[4] = (g[1]*g[6]-g[3]^2)/det     -- g^22
      gContraOut[5] = (g[2]*g[3]-g[1]*g[5])/det  -- g^23 = g^32
      gContraOut[6] = (g[1]*g[4]-g[2]^2)/det     -- g^33
   end
end

-- internal function to copy physical coordinates into a vector
-- instead of returning it
function MappedCart:_mapc2p_vec(xc, xp)
   if self._rdim == 1 then
      xp[1] = self._mapc2p(xc)
   elseif self._rdim == 2 then
      xp[1], xp[2] = self._mapc2p(xc)
   elseif self._rdim == 3 then
      xp[1], xp[2], xp[3] = self._mapc2p(xc)
   end
end

-- get cell center in physical coordinates
function MappedCart:cellCenterPhys(xp)
   self.super.cellCenter(self, self._xc)
   self:_mapc2p_vec(self._xc, xp)
end

function MappedCart:write(fName)
   -- create a grid over nodes and a field to store nodal coordinates
   local cells, lower, upper = {}, {}, {}
   for d = 1, self:ndim() do
      cells[d] = self:numCells(d)+1 -- one more layer of nodes than cells
      -- this ensures cell-center of nodal grid lie at nodes of
      -- original grid
      lower[d] = self:lower(d) - 0.5*self:dx(d)
      upper[d] = self:upper(d) + 0.5*self:dx(d)
   end
   -- WILL NEED TO MAKE THIS WORK IN PARALLEL .. EVENTUALLY
   local grid = RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
   }
   local nodalCoords = DataStruct.Field {
      onGrid = grid,
      numComponents = self._rdim,
   }

   local xnc, xnp = Lin.Vec(self:ndim()), Lin.Vec(self._rdim)

   local localRange = nodalCoords:localRange()
   local indexer = nodalCoords:genIndexer()
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)

      grid:cellCenter(xnc) -- nodal coordinate in computational space
      self:_mapc2p_vec(xnc, xnp) -- nodal coordinate in physical space

      local nitr = nodalCoords:get(indexer(idx))
      for d = 1, self._rdim do
	 nitr[d] = xnp[d]
      end
   end
   nodalCoords:write(fName)
end

return MappedCart
