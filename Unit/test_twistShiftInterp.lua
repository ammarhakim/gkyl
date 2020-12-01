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
local Lin        = require "Lib.Linalg"
local lume       = require "Lib.lume"
local root       = require "sci.root"

-- Create two fields on a 2D grid. Then interpolate one field onto the other
-- but with a shift in y that is a function of x (assuming periodicity in y).

local polyOrder    = 1
local lower        = {-2.0, -1.50}
local upper        = { 2.0,  1.50}
local numCells     = {1, 10}
local periodicDirs = {1}

local grid = Grid.RectCart {
   lower = lower,
   upper = upper,
   cells = numCells,
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
local fldSrc = DataStruct.Field {
   onGrid        = grid,
   numComponents = basis:numBasis(),
   ghost         = {1, 1},
   metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
}
local fldDest = DataStruct.Field {
   onGrid        = grid,
   numComponents = basis:numBasis(),
   ghost         = {1, 1},
   metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
}
local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn)
                 local x, y       = xn[1], xn[2]
                 local muX, muY   = 0., 0.
                 local sigX, sigY = 0.5, 0.3
                 return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
              end
}
project:advance(0., {}, {fldSrc})
--fldSrc:write("fldSrc.bp")

-- Create a 1D grid and project the function that determines the shift.
-- In flux-tube gyrokinetics this shift is a function of the magnetic
-- safety profile, something like yShift = L_z*C_y(x)*q(x).
local grid1D = Grid.RectCart {
   lower = {lower[1]},
   upper = {upper[1]},
   cells = {numCells[1]},
   periodicDirs = {},
}
local basis1D = Basis.CartModalSerendipity { ndim = grid1D:ndim(), polyOrder = polyOrder }
local yShift = DataStruct.Field {
   onGrid        = grid1D,
   numComponents = basis1D:numBasis(),
   ghost         = {1, 1},
   metaData      = {polyOrder = basis1D:polyOrder(), basisType = basis1D:id()},
}
local project1D = Updater.EvalOnNodes {
   onGrid   = grid1D,
   basis    = basis1D,
   evaluate = function(t, xn)
                 local x  = xn[1]
                 local s0 = 0.25
                 return 1./(1.+s0*x)
              end
}
project1D:advance(0., {}, {yShift})
yShift:write("yShift.bp")

-- Determine the cells that each cell in fldDest needs.
local numBasis1D = basis1D:numBasis()
local basis1Dev  = Lin.Vec(numBasis1D)

local idx1D      = Lin.IntVec(1)
local indexer    = fldSrc:genIndexer()
local indexer1D  = yShift:genIndexer()
local localRange = fldSrc:localRange()
local srcPtr, destPtr, yShiftPtr = fldSrc:get(1), fldDest:get(1), yShift:get(1)
local xc, dx     = Lin.Vec(2), Lin.Vec(2)
local idxP, dxP  = Lin.IntVec(2), Lin.Vec(2)
local stepFrac   = {0.1,0.1}          -- Size of step taken around the boundary, as a fraction of cell length.
local srcCells   = {}
local numSteps   = {math.ceil(1/stepFrac[1]),math.ceil(1/stepFrac[2])}
local wrapNum = function (val, lower, upper, pickUpper)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   local pickUpper = pickUpper==nil and false or pickUpper
   local L         = upper - lower
   local disp      = (val - lower) % L
   if math.abs(disp - GKYL_MIN_DOUBLE) > 0. then
      return lower + (L + disp) % L
   else
      local mid = 0.5*(lower+upper)
      if pickUpper and val>mid then
         return upper
      else
         return lower + (L + disp) % L
      end
   end
end
for idx in localRange:rowMajorIter() do

   -- One way of determining which cells this cell is magnetically connected to
   -- is to take points around the boundary of the cell and calculate the shift from those
   -- points, then find out which cells own the shifted points. This is perhaps limited to
   -- grid that are not refined in the z-direction; in that case we would have to consider
   -- more than just cell-boundary points.
   -- For now simply take small steps around the boundaries of the 2D cell, compute
   -- the shifted points and their owner at each step.

   grid:setIndex(idx)
   grid:cellCenter(xc)
   grid:getDx(dx)

   local cellLower  = {grid:cellLowerInDir(1),grid:cellLowerInDir(2)}
   local stepSize   = {grid:dx(1)/numSteps[1],grid:dx(2)/numSteps[2]}

   idx:copyInto(idxP)

   idx1D[1] = idx[1]
   yShift:fill(indexer1D(idx1D), yShiftPtr)
   local srcCells = {}
   for dC = 1,2 do            -- Loop over x=const and y=const boundaries.
      for xS = 1,2 do         -- Loop over lower/upper boundary.

         -- Do the first point separately (it requires pickLower=false in findCell).
         local evPoint = {cellLower[1], cellLower[2]}
         evPoint[dC]   = evPoint[dC]+(xS-1)*grid:dx(dC) 

         basis1D:evalBasis({(evPoint[1]-grid:cellCenterInDir(1))/(0.5*grid:dx(1))}, basis1Dev)
         local yShiftB = 0.   -- y-shift at boundary.
         for k = 1, numBasis1D do
            yShiftB = yShiftB + yShiftPtr[k]*basis1Dev[k]
         end

         local idxShifted0 = {nil,nil}
         grid:findCell({evPoint[1],wrapNum(evPoint[2]+yShiftB,grid:lower(2),grid:upper(2),true)},
                       idxShifted0,false,{idx[1],nil})
         if lume.findTable(srcCells,idxShifted0)==nil then table.insert(srcCells,idxShifted0) end
         
         local stepDim, newP = dC % 2+1, {nil,nil}
         newP[dC] = evPoint[dC]
         for sI = 1, numSteps[stepDim] do

            newP[stepDim] = evPoint[stepDim] + sI*stepSize[stepDim]
         
            -- Evaluate yShift at this point.
            basis1D:evalBasis({(newP[1]-grid:cellCenterInDir(1))/(0.5*grid:dx(1))}, basis1Dev)
            local yShiftB = 0.   -- y-shift at boundary.
            for k = 1, numBasis1D do
               yShiftB = yShiftB + yShiftPtr[k]*basis1Dev[k]
            end

            -- Find the index of the cell that owns the shifted point.
            local idxShifted = {nil,nil}
            grid:findCell({newP[1],wrapNum(newP[2] + yShiftB,grid:lower(2),grid:upper(2),true)},
                           idxShifted,true,{idx[1],nil})
            if lume.findTable(srcCells,idxShifted)==nil then table.insert(srcCells,idxShifted) end
         end
      end
   end

end
