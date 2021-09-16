-- Gkyl ------------------------------------------------------------------------
--
-- Test the interpolation needed for twist-shift BCs in gyrokinetics.
--
-- Create two fields on a 2D grid, the donor field and the target field. Then
-- the donor field with a shift in y (that may be a function of x) assuming
-- periodicity in y to obtain the target field.
-- 
-- This test uses a Gaussian donor field and a linear shift.
--
-- Some pgkyl commands to check the results:
--   pgkyl rt-twistShift-2x_fldDo.bp -l '$f(x,y)$' rt-twistShift-2x_fldTar.bp -l '$g(x,y)=f(x,y-\mathcal{S}(x))$' rt-twistShift-2x_fldDoBack.bp -l '$g(x,y+\mathcal{S}(x))$' pl
--   pgkyl rt-twistShift-2x_fldDo.bp -l '$f(x,y)$' rt-twistShift-2x_fldTar.bp -l '$g(x,y)=f(x,y-\mathcal{S}(x))$' rt-twistShift-2x_fldDoBack.bp -l '$g(x,y+\mathcal{S}(x))$' interp pl -b
--   pgkyl rt-twistShift-2x_fldDo_intV.bp rt-twistShift-2x_fldTar_intV.bp ev 'f[0] f[1] -' pr
--   pgkyl rt-twistShift-2x_fldDo_intV.bp rt-twistShift-2x_fldDoBack_intV.bp ev 'f[0] f[1] -' pr
-- 
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local polyOrder       = 1
local lower           = {-2.0, -1.50}
local upper           = { 2.0,  1.50}
local cells           = {{1,10}, {10,10}, {10,5}, {20, 10}, {40,20}, {80,40}, {160,80}, {320,160}}
local periodicDirs    = {2}
local yShiftPolyOrder = polyOrder  -- Poly order for discretizing the shift function, = polyOrder, or 1 or 2 for polyOrder=1.

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
                      return 0.6*x+1.8
                   end
-- Donor field function.
local fldDoFunc = function(t, xn)
   local x, y       = xn[1], xn[2]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.6, 0.2
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
end

-- ....................... END OF USER INPUTS (maybe) ........................... --

local yShiftBackFunc = function(t, xn) return -yShiftFunc(t, xn) end

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
   local x, y = xn[1], xn[2]
   local yS   = wrapNum(y-yShiftFunc(0,xn),{lo=lower[2],up=upper[2]},true)
   return fldDoFunc(t, {x, yS})
end

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

local basis = Basis.CartModalSerendipity { ndim = #lower, polyOrder = polyOrder }

local grids = {}
local fldDos, fldDoShifteds, fldTars = {}, {}, {}
local projects, intQuants = {}, {}
local intFldDos, intFldTars, intFldDoBacks = {}, {}, {}
local twistShiftUpds, twistShiftBackUpds = {}, {}

for gI, numCells in ipairs(cells) do 

   grids[gI] = Grid.RectCart {
      lower = lower,  cells        = numCells,
      upper = upper,  periodicDirs = periodicDirs,
   }
   local grid = grids[gI]

   local function fileName(varNm)
      return string.format(varNm .. "_Nx%dNy%dP%d_yShP%d.bp",grid:numCells(1),grid:numCells(2),polyOrder,yShiftPolyOrder)
   end
   
   fldDos[gI]        = createField(grid, basis)
   fldDoShifteds[gI] = createField(grid, basis)
   fldTars[gI]       = createField(grid, basis)
   local fldDo, fldDoShifted, fldTar = fldDos[gI], fldDoShifteds[gI], fldTars[gI]
   
   projects[gI] = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1. end,  -- Set later.
   }
   local project = projects[gI]
   -- Project donor field function onto basis.
   project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
   project:advance(0., {}, {fldDo})
   fldDo:write(fileName("fldDo"))
   -- Project shifted donor field function onto basis.
   project:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
   project:advance(0., {}, {fldDoShifted})
   fldDoShifted:write(fileName("fldDoShifted"))
   
   -- Compute the integral of the donor field.
   intQuants[gI] = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V",
   }
   local intQuant = intQuants[gI]
   intFldDos[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldDo = intFldDos[gI]
   intQuant:advance(0., {fldDo}, {intFldDo})
   intFldDo:write(fileName("fldDo_intV"), 0., 0)
   
   twistShiftUpds[gI] = Updater.TwistShiftBC {
      onGrid = grid,   yShiftFunc      = yShiftFunc, 
      basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
   }
   local twistShiftUpd = twistShiftUpds[gI]
   
   -- Apply the forward shift.
   local t1 = os.clock()
   twistShiftUpd:_advance2x(0., {fldDo}, {fldTar})
   local t2 = os.clock()
   --io.write("Forward shift time: ", t2-t1, " s\n")
   
   fldTar:write(fileName("fldTar"))
   
   -- Compute the integral of the target field.
   intFldTars[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldTar = intFldTars[gI] 
   intQuant:advance(0., {fldTar}, {intFldTar})
   intFldTar:write(fileName("fldTar_intV"), 0., 0)
   
   
   -- ............... SHIFT BACK .................. --
   twistShiftBackUpds[gI] = Updater.TwistShiftBC {
      onGrid = grid,   yShiftFunc      = yShiftBackFunc, 
      basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
   }
   local twistShiftBackUpd = twistShiftBackUpds[gI]
   
   fldDo:clear(0.)
   local t1 = os.clock()
   -- Apply the backward shift.
   twistShiftBackUpd:_advance2x(0., {fldTar}, {fldDo})
   local t2 = os.clock()
   --io.write("Backward shift time: ", t2-t1, " s\n")
   
   fldDo:write(fileName("fldDoBack"))
   
   -- Compute the integral of the new target field.
   intFldDoBacks[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldDoBack = intFldDoBacks[gI]
   intQuant:advance(0., {fldDo}, {intFldDoBack})
   intFldDoBack:write(fileName("fldDoBack_intV"), 0., 0)

end
