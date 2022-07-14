-- Gkyl ------------------------------------------------------------------------
--
-- Mapped Cartesian grids: computational space is rectangular but
-- physical space need not be.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local DataStruct = require "DataStruct"
local Lin        = require "Lib.Linalg"
local Proto      = require "Lib.Proto"
local RectCart   = require "Grid.RectCart"
local diff       = require "sci.diff-recursive"
local diff1      = require "sci.diff"
local math       = require("sci.math").generic

-- MappedCartGrid --------------------------------------------------------------
--
-- Grid determined by a coordinate mapping between uniform rectangular
-- computational space to physical space.
--------------------------------------------------------------------------------

local MappedCart = Proto(RectCart) -- Extends RectCart.
-- ctor

function MappedCart:init(tbl)
   MappedCart.super.init(self, tbl)

   -- Function that maps the computational coordinates to
   -- cartesian coordinates of the form 
   -- X, Y, Z = f({xcomp1, xcomp2, xcomp3})
   self._mapc2p = tbl.mapc2p

   -- Domain of the higher dimensional world this simulation lives in.
   -- These optional parameters are used to call mapc2p with more computational
   -- coordinates than the simulation actually has (to compute metric).
   self._useWorld   = false
   self._worldCoord = nil
   if tbl.world then
      self._worldDim   = tbl.world["dim"]
      local extraCoord = tbl.world["evaluateAt"]
      assert(self._worldDim > self:ndim(), "MappedCart: 'dim' in 'world' must be larger than dimensionality of simulation.")

      self._addDim = 0
      for nm, _ in pairs(extraCoord) do self._addDim=self._addDim+1 end
      assert(self._addDim < self._worldDim,
             "MappedCart: evaluateAt entries must be fewer than dim.")
      self._useWorld = true

      -- Create a way to complement computational coordinates
      -- with additional coordinates
      self._compIdx  = {}
      self._addIdx   = {}
      self._addCoord = {}

      local cLabels, cTrans = {"x","y","z","vx","vy","vz"}, {}
      for d = 1, self._worldDim do cTrans[cLabels[d]]=d end
      local xcI, addI = 0, 0
      lume.setOrder(cTrans)
      for nm, v in lume.orderedIter(cTrans) do -- Order in this loop matters.
         if extraCoord[nm] then 
            addI = addI+1
            self._addIdx[addI]   = v 
            self._addCoord[addI] = extraCoord[nm]
         else
            xcI = xcI+1
            self._compIdx[xcI] = v
         end
      end

      self._worldCoord = {}
      for d = 1, self:ndim() do self._worldCoord[self._compIdx[d]] = self._lower[d] end
      for d = 1, self._addDim do self._worldCoord[self._addIdx[d]] = self._addCoord[d] end
      self._inDim = self._worldDim                        -- Dimensionality of the input to mapc2p.
      self._rdim  = #{ self._mapc2p(self._worldCoord) }   -- Number of values mapc2p returns.
   else
      self._inDim = self:ndim()                      -- Dimensionality of the input to mapc2p.
      self._rdim  = #{ self._mapc2p(self._lower) }   -- Number of values mapc2p returns.
   end

   self._xc = Lin.Vec(self._inDim)
   self._d1, self._d2, self._d3 = Lin.Vec(self._inDim), Lin.Vec(self._inDim), Lin.Vec(self._inDim)
end

--function MappedCart:id() return "mapped" end
function MappedCart:rdim() return self._rdim end

function MappedCart:mapc2p(xc)
   return self._mapc2p(xc)
end

-- Internal methods, not to be used directly by user.
function MappedCart:_calcMetric_1d_r1(xc, g)
   local d1 = self._d1
   d1[1] = diff.derivt(self._mapc2p, 1)(xc)

   g[1] = d1[1]^2
end
function MappedCart:_calcMetric_1d_r2(xc, g)
   local d1, d2 = self._d1, self._d2
   d1[1], d2[1] = diff.derivt(self._mapc2p, 1)(xc)

   g[1] = d1[1]^2 + d2[1]^2
end
function MappedCart:_calcMetric_1d_r3(xc, g)
   local d1, d2, d3 = self._d1, self._d2, self._d3
   d1[1], d2[1], d3[1] = diff.derivt(self._mapc2p, 1)(xc)

   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2
end
function MappedCart:_calcMetric_2d_r2(xc, g)
   local d1, d2 = self._d1, self._d2
    
   d1[1], d2[1] = diff.derivt(self._mapc2p, 1)(xc)
   d1[2], d2[2] = diff.derivt(self._mapc2p, 2)(xc)

   g[1] = d1[1]^2 + d2[1]^2         -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] -- g_12 = g_21
   g[3] = d1[2]^2 + d2[2]^2         -- g_22
end
function MappedCart:_calcMetric_2d_r3(xc, g)
   local d1, d2, d3 = self._d1, self._d2, self._d3
   d1[1], d2[1] = diff.derivt(self._mapc2p, 1)(xc)
   d1[2], d2[2] = diff.derivt(self._mapc2p, 2)(xc)
   d1[3], d2[3] = diff.derivt(self._mapc2p, 3)(xc)

   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2              -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] + d3[1]*d3[2]  -- g_12 = g_21
   g[3] = d1[2]^2 + d2[2]^2 + d3[2]^2              -- g_22
end
function MappedCart:_calcMetric_3d(xc, g)
   local d1, d2, d3 = self._d1, self._d2, self._d3
   d1[1], d2[1], d3[1] = diff.derivt(self._mapc2p, 1)(xc)
   d1[2], d2[2], d3[2] = diff.derivt(self._mapc2p, 2)(xc)
   d1[3], d2[3], d3[3] = diff.derivt(self._mapc2p, 3)(xc)

   g[1] = d1[1]^2 + d2[1]^2 + d3[1]^2              -- g_11
   g[2] = d1[1]*d1[2] + d2[1]*d2[2] + d3[1]*d3[2]  -- g_12 = g_21
   g[3] = d1[1]*d1[3] + d2[1]*d2[3] + d3[1]*d3[3]  -- g_13 = g_31
   g[4] = d1[2]^2 + d2[2]^2 + d3[2]^2              -- g_22
   g[5] = d1[2]*d1[3] + d2[2]*d2[3] + d3[2]*d3[3]  -- g_23 = g_32
   g[6] = d1[3]^2 + d2[3]^2 + d3[3]^2              -- g_33
end

-- Computes (covariant) metric tensor g_ij.
function MappedCart:calcMetric(xc, gOut)
   if self._useWorld then
      -- Complement the computational coordinates xc w/ dimensions not simulated.
      for d = 1, self:ndim() do self._worldCoord[self._compIdx[d]] = xc[d] end
   else
      self._worldCoord = xc
   end

   if self._inDim == 1 then
      if self._rdim == 1 then
	 self:_calcMetric_1d_r1(self._worldCoord, gOut)
      elseif self._rdim == 2 then
	 self:_calcMetric_1d_r2(self._worldCoord, gOut)
      elseif self._rdim == 3 then
	 self:_calcMetric_1d_r3(self._worldCoord, gOut)
      end
   elseif self._inDim == 2 then
      if self._rdim == 2 then
	 self:_calcMetric_2d_r2(self._worldCoord, gOut)
      else
	 self:_calcMetric_2d_r3(self._worldCoord, gOut)
      end
   elseif self._inDim == 3 then
      self:_calcMetric_3d(self._worldCoord, gOut)
   else
      assert(false, "MappedCart does not support more than 3 dimensions!")
   end
end

-- Calculates gradients for jacobian ONLY set up for xyz case.
function MappedCart:mapDiff(xv)
   -- Compute gradients.
   local grad1 = diff1.gradientf(function (xc) local x, _, _ = self:mapc2p(xc) return x end, 2)
   local dx = Lin.Vec(2)
   grad1(xv, dx)
   local grad2 = diff1.gradientf(function (xc) local _, y, _ = self:mapc2p(xc) return y end, 2)
   local dy = Lin.Vec(2)
   grad2(xv, dy)
   local grad3 = diff1.gradientf(function (xc) local _, _, z = self:mapc2p(xc) return z end, 2)
   local dz = Lin.Vec(2)
   grad3(xv, dz)

   -- Separate derivatives.
   local dxdxc = dx[1]
   local dydxc = dy[1]
   local dzdxc = dz[1]
   local dxdyc = dx[2]
   local dydyc = dy[2]
   local dzdyc = dz[2]

   return dxdxc, dydxc, dzdxc, dxdyc, dydyc, dzdyc
end

-- Computes jacobian components
function MappedCart:calcDiffLen(xcv, hOut)
    -- Get derivatives, physical coordinates.
    local dxdxc, dydxc, dzdxc, dxdyc, dydyc, dzdyc = self:mapDiff(xcv)
    -- Assign to input vector.
    hOut[1], hOut[2], hOut[3], hOut[4], hOut[5], hOut[6] = dxdxc, dydxc, dzdxc, dxdyc, dydyc, dzdyc
end

-- Compute jacobian = det(g_ij)^(1/2).
function MappedCart:calcJacobian(xc)
   local g = {}
   self:calcMetric(xc, g)
   local jacobian
   if self._inDim == 1 then
      jacobian = math.sqrt(g[1]*g[1])
   elseif self._inDim == 2 then
      jacobian = math.sqrt(g[1]*g[3] - g[2]*g[2]) 
   elseif self._inDim == 3 then
      jacobian = math.sqrt(-g[3]^2*g[4] + 2*g[2]*g[3]*g[5] - g[1]*g[5]^2 - g[2]^2*g[6] + g[1]*g[4]*g[6])
   end
   return jacobian
end  

-- Computes contravariant metric tensor g^ij = (g_ij)^-1
function MappedCart:calcContraMetric(xc, gContraOut)
   local g = {}
   self:calcMetric(xc, g)
   if self._inDim == 1 then
      gContraOut[1] = 1/g[1]
   elseif self._inDim == 2 then
      local det = self:calcJacobian(xc)^2
      gContraOut[1] = g[3]/det   -- g^11
      gContraOut[2] = -g[2]/det  -- g^12
      gContraOut[3] = g[1]/det   -- g^22
   elseif self._inDim == 3 then
      local det = self:calcJacobian(xc)^2
      gContraOut[1] = (g[4]*g[6]-g[5]^2)/det     -- g^11
      gContraOut[2] = (g[3]*g[5]-g[2]*g[6])/det  -- g^12 = g^21
      gContraOut[3] = (g[2]*g[5]-g[3]*g[4])/det  -- g^13 = g^31
      gContraOut[4] = (g[1]*g[6]-g[3]^2)/det     -- g^22
      gContraOut[5] = (g[2]*g[3]-g[1]*g[5])/det  -- g^23 = g^32
      gContraOut[6] = (g[1]*g[4]-g[2]^2)/det     -- g^33
   end
end

-- Internal function to copy physical coordinates into
-- a vector instead of returning separate outputs.
function MappedCart:_mapc2p_vec(xc, xp)
   if self._rdim == 1 then
      xp[1] = self._mapc2p(xc)
   elseif self._rdim == 2 then
      xp[1], xp[2] = self._mapc2p(xc)
   elseif self._rdim == 3 then
      xp[1], xp[2], xp[3] = self._mapc2p(xc)
   end
end

-- Get cell center in physical coordinates.
function MappedCart:cellCenterPhys(xp)
   self.super.cellCenter(self, self._xc)
   self:_mapc2p_vec(self._xc, xp)
end

function MappedCart:write(fName)
   -- Write a file containing the grid node coordinates.

   -- Create a grid over nodes and a field to store nodal coordinates.
   local cells, lower, upper = {}, {}, {}
   for d = 1, self:ndim() do
      cells[d] = self:numCells(d)+1   -- One more layer of nodes than cells.
      -- This ensures cell-center of nodal grid lie at nodes of original grid
      lower[d] = self:lower(d) - 0.5*self:dx(d)
      upper[d] = self:upper(d) + 0.5*self:dx(d)
   end
   -- Create a grid of nodes.
   local grid = RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      decomposition = self.decomp,
   }
   local nodalCoords = DataStruct.Field {
      onGrid        = grid,
      numComponents = self._rdim,
   }

   local xnc, xnp = Lin.Vec(self._inDim), Lin.Vec(self._rdim)

   local localRange = nodalCoords:localRange()
   local indexer    = nodalCoords:genIndexer()
   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)

      grid:cellCenter(xnc)         -- Nodal coordinate in computational space.

      if self._useWorld then
         -- Complement the computational coordinates xc w/ dimensions not simulated.
         for d = 1, self:ndim() do self._worldCoord[self._compIdx[d]] = xnc[d] end
      else
         self._worldCoord = xnc
      end

      self:_mapc2p_vec(self._worldCoord, xnp)   -- Nodal coordinate in physical space.

      local nPtr = nodalCoords:get(indexer(idx))
      for d = 1, self._rdim do nPtr[d] = xnp[d] end
   end
   nodalCoords:write(fName)
end

return MappedCart
