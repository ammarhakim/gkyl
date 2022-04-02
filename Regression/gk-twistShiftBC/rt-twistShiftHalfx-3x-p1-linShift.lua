-- Gkyl ------------------------------------------------------------------------
--
-- Test the interpolation needed for twist-shift BCs in gyrokinetics over
-- half the x domain only.
--
-- This tests the 3x twist shift by taking a field that is only non-zero in
-- the inner cells, and using the twistShift update to populate its z ghost
-- cells.
--
-- In order to check that the integral of that ghost cell is the same as that
-- of the skin cell it originated from, we need to create a field that projects
-- the donor function only in the ghost cell or only in the skin cell.
-- 
-- This test uses a Gaussian donor field and a linear shift.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Range      = require "Lib.Range"

local polyOrder       = 1
local lower           = {-2.0, -1.50, -3.0}
local upper           = { 2.0,  1.50,  3.0}
local cells           = {{1,10,4}, {10,10,4}, {20,20,4}, {40,40,4}, {80,80,4}}  -- Code below assumes Nz is the same for all.
local periodicDirs    = {2}
local yShiftPolyOrder = polyOrder  -- Polyorder for discretizing the shift function, = polyOrder, or 1 or 2 for polyOrder=1.

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
                      return 0.3*x+1.4
                   end
-- Donor field function.
local fldDoFunc = function(t, xn)
   local x, y       = xn[1], xn[2]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
end

-- ....................... END OF USER INPUTS (maybe) ........................... --

local dz = (upper[3]-lower[3])/cells[1][3]

local shiftFuncs = {lower=yShiftFunc, upper=function(t,xn) return -yShiftFunc(t,xn) end}

-- Donor field function that only allows projection in the z skin cell.
local fldDoFuncZskinOnly = {
   lower = function(t, xn)
      local x, y, z = xn[1], xn[2], xn[3]
      if ((   x > lower[1]) and (   x < upper[1])) and
         ((   y > lower[2]) and (   y < upper[2])) and
         ((   z > lower[3]) and (   z < lower[3]+dz)) then
         return fldDoFunc(t, xn)
      else
         return 0.
      end
   end,
   upper = function(t, xn)
      local x, y, z = xn[1], xn[2], xn[3]
      if ((   x > lower[1])    and (   x < upper[1])) and
         ((   y > lower[2])    and (   y < upper[2])) and
         ((   z > upper[3]-dz) and (   z < upper[3])) then
         return fldDoFunc(t, xn)
      else
         return 0.
      end
   end
}

local wrapNum = function (val, lims, pickUpper)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   local lower, upper = lims.lo, lims.up
   local L        = upper - lower
   local disp     = (val - lower) % L
   local newCoord = lower + (L + disp) % L
   local eps      = 1.e-12
   if ( (lower-eps < newCoord and newCoord < lower + eps) or
        (upper-eps < newCoord and newCoord < upper + eps) ) then
      if pickUpper then 
         return upper
      else
         return lower
      end
   else
      return newCoord
   end
end

-- Shifted donor field function (should be the same as fldDoFunc but shifted in y).
local fldDoShiftedFunc = function(t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local yS = wrapNum(y-yShiftFunc(0,xn),{lo=lower[2],up=upper[2]},true)
   return fldDoFunc(t, {x, yS, z})
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
   }
   fld:clear(0.)
   return fld
end

function getGhostRange(edge, dir, global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if edge == "lower" then
      uv[dir] = global:lower(dir)-1   -- For ghost cells on "left".
   else
      lv[dir] = global:upper(dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

local basis = Basis.CartModalSerendipity { ndim = #lower, polyOrder = polyOrder }

local edges     = {"lower","upper"}
local dualEdges = {"upper","lower"}

local grids = {}
local fldDos, fldDoShifteds, fldTars = {}, {}, {}
local fldDoNoGhosts = {}
local projects, projectOnGhosts, intQuants = {}, {}, {}
local intFldDos, intFldTars, intFldDoBacks = {}, {}, {}
local twistShiftUpds, twistShiftBackUpds = {}, {}

for gI, numCells in ipairs(cells) do 

   grids[gI] = Grid.RectCart {
      lower = lower,  cells        = numCells,
      upper = upper,  periodicDirs = periodicDirs,
   }
   local grid = grids[gI]

   local function fileName(varNm)
      return string.format(varNm .. "_Nx%dNy%dNz%dP%d_yShP%d.bp",grid:numCells(1),grid:numCells(2),grid:numCells(3),polyOrder,yShiftPolyOrder)
   end
   
   fldDos[gI]        = createField(grid, basis)
   fldDoShifteds[gI] = createField(grid, basis)
   fldTars[gI]       = createField(grid, basis)
   fldDoNoGhosts[gI] = createField(grid, basis)
   local fldDo, fldDoShifted, fldTar = fldDos[gI], fldDoShifteds[gI], fldTars[gI]
   local fldDoNoGhost = fldDoNoGhosts[gI]
   
   projects[gI] = Updater.ProjectOnBasis {
      onGrid = grid,   evaluate = function(t, xn) return 1. end,
      basis  = basis,
   }
   local project = projects[gI]
   projectOnGhosts[gI] = Updater.ProjectOnBasis {
      onGrid = grid,   evaluate = function(t, xn) return 1. end,
      basis  = basis,  onGhosts = true,
   }
   local projectOnGhost = projectOnGhosts[gI]

   -- Project donor field function onto basis.
   projectOnGhost:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
   projectOnGhost:advance(0., {}, {fldDo})
   fldDo:write(fileName("fldDo"))
   -- Project shifted donor field function onto basis.
   projectOnGhost:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
   projectOnGhost:advance(0., {}, {fldDoShifted})
   fldDoShifted:write(fileName("fldDoShifted"))

   intFldDos[gI]  = {}
   intFldTars[gI] = {}
   local intFldDo  = intFldDos[gI]
   local intFldTar = intFldTars[gI] 

   intQuants[gI]  = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V",
   }
   local intQuant = intQuants[gI]

   twistShiftUpds[gI] = {}
   local twistShiftUpd = twistShiftUpds[gI]

   for i, edge in ipairs(edges) do

      -- Create ghost range with half the x-domain over which to apply BCs:
      local global, globalExt = fldTar:globalRange(), fldTar:globalExtRange()
      local localExtRange = fldTar:localExtRange()
      local ghostRangeAll = localExtRange:intersect(getGhostRange(edge, 3, global, globalExt))
      local ghostRange    = ghostRangeAll:shorten(1,1+numCells[1]/2) -- Lower x half.
      --local ghostRange    = ghostRangeAll:shortenFromBelow(1,1+numCells[1]/2) -- Upper x half.

      -- Project the donor field function into the target field but not
      -- in the ghost cells. Twist-shift will fill the ghost cells.
      project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
      project:advance(0., {}, {fldDoNoGhost})
      fldTar:copy(fldDoNoGhost)
   
      local dualEdge = dualEdges[i]
      -- Compute the velocity moments of the donor field in the upper ghost cells only,
      -- to be compared with the moments of the target field in the lower ghost cells.
      projectOnGhost:setFunc(function(t,xn) return fldDoFuncZskinOnly[dualEdge](t,xn) end)
      projectOnGhost:advance(0., {}, {fldDo})
      fldDo:write(fileName("fldDo_"..dualEdge), 0, 0)
   
      -- Compute the integral of the donor field.
      intFldDo[edge] = DataStruct.DynVector { numComponents = 1, }
      intQuant:advance(0., {fldDo}, {intFldDo[edge]})
      intFldDo[edge]:write(fileName("fldDo_intV_"..dualEdge), 0., 0)
   
      twistShiftUpd[edge] = Updater.TwistShiftBC {
         onGrid = grid,   yShiftFunc      = shiftFuncs[edge], 
         basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
         edge   = edge,   ghostRange      = ghostRange,
      }
   
      -- Apply the forward shift.
      local t1 = os.clock()
      twistShiftUpd[edge]:_advance3xInPlace(0., {}, {fldTar})
      local t2 = os.clock()
      --io.write([edge].." shift time: ", t2-t1, " s\n")

      -- Subtract the donor field in the inner cells so we only compute
      -- ghost cell moments.
      fldTar:accumulate(-1., fldDoNoGhost)
      fldTar:write(fileName("fldTar_"..edge), 0., 0, true)
   
      -- Compute the integral of the target field.
      intFldTar[edge] = DataStruct.DynVector { numComponents = 1, }
      intQuant:advance(0., {fldTar}, {intFldTar[edge]})
      intFldTar[edge]:write(fileName("fldTar_intV_"..edge), 0., 0)
   
   end
end
