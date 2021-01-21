-- Gkyl ------------------------------------------------------------------------
--
-- Test the interpolation needed for twist-shift BCs in gyrokinetics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
-- The following are needed for projecting onto the basis.
local SerendipityNodes     = require "Lib.SerendipityNodes"
GKYL_EMBED_INP             = false

-- Create two fields on a 2D grid. Then interpolate one field onto the other
-- but with a shift in y that is a function of x (assuming periodicity in y).

local polyOrder       = 1
local lower           = {-2.0, -1.50, -3.0}
local upper           = { 2.0,  1.50,  3.0}
local numCells        = {10, 10, 4}
local periodicDirs    = {2}
local yShiftPolyOrder = 1

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
   }
   return fld
end

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

local grid = Grid.RectCart {
   lower        = lower,
   upper        = upper,
   cells        = numCells,
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
local fldDo        = createField(grid, basis)
local fldDoShifted = createField(grid, basis)
local fldTar       = createField(grid, basis)

-- Create a 1D grid and project the function that determines the shift.
-- In flux-tube gyrokinetics this shift is a function of the magnetic
-- safety profile, something like yShift = L_z*C_y(x)*q(x).
local polyOrder1D = polyOrder
local grid1D = Grid.RectCart {
   lower        = {lower[1]},
   upper        = {upper[1]},
   cells        = {numCells[1]},
   periodicDirs = {},
}
local basis1D = Basis.CartModalSerendipity { ndim = grid1D:ndim(), polyOrder = polyOrder1D }

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
--                      return 1./(1.+0.25*x)
                      return -0.3*x+0.97
--                      return 0.1*(x+2.)^2+0.2*x+0.6
                   end

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn) return 1. end,
   projectOnGhosts = true,
}
-- Donor field.
local fldDoFunc = function(t, xn)
   local x, y, z    = xn[1], xn[2], xn[3]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
--   return math.exp(-((y-muY)^2)/(2.*(sigY^2)))
--   return math.sin((2.*math.pi/3.)*y)
--   if y < 0. then
--      return 0.
--   else
--      return 1.
--   end
end
-- Shifted donor field.
local fldDoShiftedFunc = function(t, xn)
   local x, y, z    = xn[1], xn[2], xn[3]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((wrapNum(y+yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true)-muY)^2)/(2.*(sigY^2)))
--   return math.exp(-((wrapNum(y+yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true)-muY)^2)/(2.*(sigY^2)))
--   return math.sin((2.*math.pi/3.)*((wrapNum(y+yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true))))
--   if wrapNum(y+yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true) < 0. then
--      return 0.
--   else
--      return 1.
--   end
end

-- Project donor field function onto basis.
project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
project:advance(0., {}, {fldDo})
fldDo:write("fldDo.bp", 0., 0, true)
-- Project shifted donor field function onto basis.
project:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
project:advance(0., {}, {fldDoShifted})
fldDoShifted:write("fldDoShifted.bp", 0., 0, true)

-- Project the function in the target field but not in the ghost
-- cells. Twist-shift will fill the ghost cells.
local projectNoGhosts = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn) return 1. end,
}
projectNoGhosts:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
projectNoGhosts:advance(0., {}, {fldTar})

local intQuant = Updater.CartFieldIntegratedQuantCalc {
   onGrid        = grid,
   basis         = basis,
   numComponents = 1,
   quantity      = "V",
}
local intFldDo = DataStruct.DynVector { numComponents = 1, }
intQuant:advance(0., {fldDo}, {intFldDo})
intFldDo:write("intFldDo.bp",0., 0)

local twistShiftUpd = Updater.TwistShift {
   onGrid          = grid,
   basis           = basis, 
   yShiftFunc      = yShiftFunc, 
   yShiftPolyOrder = yShiftPolyOrder,
   edge            = "lower",
}

local t1 = os.clock()
twistShiftUpd:_advance(0., {}, {fldTar})
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")

fldTar:write("fldTar.bp", 0., 0, true)

local intFldTar = DataStruct.DynVector { numComponents = 1, }
intQuant:advance(0., {fldTar}, {intFldTar})
intFldTar:write("intFldTar.bp", 0., 0)
