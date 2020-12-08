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
-- The following are needed for projecting onto the basis.
local SerendipityNodes = require "Lib.SerendipityNodes"
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local TwistShiftInterpDecl = require "Updater.interpolateCalcData.TwistShiftInterpModDecl"
GKYL_EMBED_INP = false

-- Create two fields on a 2D grid. Then interpolate one field onto the other
-- but with a shift in y that is a function of x (assuming periodicity in y).

local polyOrder    = 1
local lower        = {-2.0, -1.50}
local upper        = { 2.0,  1.50}
local numCells     = {1, 10}
local periodicDirs = {2}

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

-- Create a 1D grid and project the function that determines the shift.
-- In flux-tube gyrokinetics this shift is a function of the magnetic
-- safety profile, something like yShift = L_z*C_y(x)*q(x).
local polyOrder1D = polyOrder
local grid1D = Grid.RectCart {
   lower = {lower[1]},
   upper = {upper[1]},
   cells = {numCells[1]},
   periodicDirs = {},
}
local basis1D = Basis.CartModalSerendipity { ndim = grid1D:ndim(), polyOrder = polyOrder1D }
local yShift = DataStruct.Field {
   onGrid        = grid1D,
   numComponents = basis1D:numBasis(),
   ghost         = {1, 1},
   metaData      = {polyOrder = basis1D:polyOrder(), basisType = basis1D:id()},
}
local project1D = Updater.EvalOnNodes {
   onGrid   = grid1D,
   basis    = basis1D,
   evaluate = function(t, xn) return 1. end   -- Set below.
}
local yShiftFunc = function(x)
--                      local s0 = 0.25
--                      return 1./(1.+s0*x)
                      return -0.3*x+1.
                   end
project1D:setFunc(function(t, xn) return yShiftFunc(xn[1]) end)
project1D:advance(0., {}, {yShift})
yShift:write("yShift.bp")

-- Invert yShift(x).
local yShiftInvFunc = function(ySIn)
   local function lossF(xIn)
      return ySIn - yShiftFunc(xIn)
   end
   local function rootStop(eps)   -- Root-finding stopping criteria.
      return function(x, y, xl, xu, yl, yu)
         if math.abs(y) < eps then return true else return false end
      end
   end
   local tol   = 1.e-13
   return root.ridders(lossF, lower[1], upper[1], rootStop(tol))
end
local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn)
                 local x, y       = xn[1], xn[2]
--                 local muX, muY   = 0., 0.
--                 local sigX, sigY = 0.5, 0.3
--                 return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
--                 return 1.
--                 return yShiftFunc(x)*math.sin((2.*math.pi/(upper[2]-lower[2]))*y)
                  return y-yShiftFunc(x) --+2  --yShiftInvFunc(yShiftFunc(x))
              end
}
project:advance(0., {}, {fldSrc})
fldSrc:write("fldSrc.bp")
local intQuant = Updater.CartFieldIntegratedQuantCalc {
   onGrid        = grid,
   basis         = basis,
   numComponents = 1,
   quantity      = "V",
}
local intFldSrc = DataStruct.DynVector {
   numComponents = 1,
}
intQuant:advance(0., {fldSrc}, {intFldSrc})
intFldSrc:write("intFldSrc.bp",0.,0)

local wrapNum = function (val, lower, upper, pickUpper)
--local wrapNum = function (val, lower, upper, idx, idxLim)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   --local pickUpper = pickUpper==nil and false or pickUpper
   local L        = upper - lower
   local disp     = (val - lower) % L
   local eps      = 1.e-10
   local newCoord = lower + (L + disp) % L
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
--   if math.abs(disp) > GKYL_MIN_DOUBLE then
--      return newCoord
--   else
--      local mid = 0.5*(lower+upper)
--      if pickUpper and val>mid then
--         return upper
--      else
--         return newCoord
--      end
--   end
end
--print(wrapNum(1.5,grid:lower(2),grid:upper(2),true,true))
--assert(false,"stop here")

local function rootStop(eps)   -- Root-finding stopping criteria.
   return function(x, y, xl, xu, yl, yu)
      if math.abs(y) < eps then return true else return false end
   end
end

-- Given a y-coordinate of the destination cell (yDest), and a y-coordinate of the
-- source cell (ySrc), find the x-coordinate of the point where the y+yShift(x) and
-- y=ySrc lines intersect. Need to provide a window in x where to search ([xBounds[1],xBounds[2]]).
local yShiftedyIntersect = function(xBounds, yDest, ySrc)
   local function lossF(xIn)
      -- Rather than wrapping yDest+yShift, it might be more robust
      -- to look for the intersection in an extended space.
      local shiftedY = yDest+yShiftFunc(xIn)
      -- For monotonically decreasing yShift.
      if ySrc <= yDest then
         return shiftedY - (grid:upper(2) + ySrc-grid:lower(2))
      else
         return shiftedY - ySrc 
      end
   end
   local tol  = 1.e-13
   -- Check if the root is bracketed by xBounds (otherwise root.ridders throws an error).
   -- We however use the lack of a root to identify the topology of the source cell region.
   local lossFl, lossFu = lossF(xBounds[1]), lossF(xBounds[2])
   if lossFl*lossFu < 0. then
      return root.ridders(lossF, xBounds[1], xBounds[2], rootStop(tol))
   else
      return nil
   end
end


-- We want to create a kernel that computes
--   int_{yL}^{yU} int_{qInvL(y)}^{qInvU(y)} f dx dy
-- where
--   yL,yU: normalized lower/upper boundary of the source cell (e.g. -1,1).
--   qInv(y): inverse of the normalized yShift.
-- For the latter consider take q(x) to be a function in x\in[-1,1] which
-- returns the yShift line in the normalized [-1,1] space.
-- Therefore, given the yShiftFunc, we need:
--   1. yShiftNormFunc(xN), which takes in the normalized x (which is defined in [-1,1])
--      computes the physical x, then computes the yShift, and returns the normalized
--      yShiftN coordinate in [-1,1] space.
--   2. yShiftNormInvFunc(yN), the inverse of yShiftNormFunc(xN). So it takes
--      a normalized y-coordinate (in [-1,1]) and returns a normalized x-coordinate.

-- Given a normalized x, the source and destination cell centers and the
-- destination cell length, compute the normalized shifted y (i.e. in [-1,1]).
local yShiftNormFunc = function(xN, y0, xcSrc, xcDest, dxDest, pickUpper)
   local xPhys = xcDest[1] + 0.5*dxDest[1]*xN
   local yS    = y0+yShiftFunc(xPhys)
--   print("xN= ",xN, "xPhys =",xPhys, "y0 =",y0, "yS =",yS)
   yS = wrapNum(yS,grid:lower(2),grid:upper(2),pickUpper)
   local ySN   = (yS-xcSrc[2])/(0.5*dxDest[2])
--   print("xN= ",xN, "xPhys =",xPhys, "y0 =",y0, "yS =",yS, "ySN =",ySN)
   return ySN
end

-- Given a normalized yShift, the source and destination cell centers and the
-- destination cell length, compute the normalzied x (i.e. in [-1,1]).
local yShiftNormInvFunc = function(ySNIn, y0In, xcSrcIn, xcDestIn, dxDestIn)
   local function lossF(xNIn)
      return ySNIn - yShiftNormFunc(xNIn, y0In, xcSrcIn, xcDestIn, dxDestIn)
   end
   local tol   = 1.e-13
   local xOutN = root.ridders(lossF, -1., 1., rootStop(tol))
   return xOutN
end
local yShiftNormInvFuncEx = function(ySNIn, y0In, xcSrcIn, xcDestIn, dxDestIn)
   local function lossF(xNIn)
      return ySNIn - yShiftNormFunc(xNIn, y0In, xcSrcIn, xcDestIn, dxDestIn)
   end
   local tol = 1.e-13
   -- Have to consider extended logical space even though the part of the y+yShift curve
   -- that is not in this cell does not get used in the integral because the fixed
   -- y limit cuts if off. But if we didn't consider extended logical space the shape of the
   -- x-limit function would not come out right.
--   print("ySN= ",ySNIn, " | y0= ",y0In)
   local xOutN = root.ridders(lossF, -2., 2., rootStop(tol))
   return xOutN
end
local yShiftNormInvFuncPartialy = function(ySNIn, y0In, xcSrcIn, xcDestIn, dxDestIn, xLim, pickUpper)
   -- y+yShift intersects the x-boundary instead of the y-boundary. This means that the y+yShift
   -- line that the defines the x-intergral limit is contained within a subset of the y-range
   -- of the source cell.
   -- This function inverts that y+yShifted(x) line. This function gets passed to a projection
   -- operation, so that the resulting DG expansion is defined in a subcell of the source cell.
   local function lossF(xNIn)
--      print("ySIn = ",ySNIn)
      return ySNIn - yShiftNormFunc(xNIn, y0In, xcSrcIn, xcDestIn, dxDestIn, pickUpper)
   end
   local tol = 1.e-12
   -- Have to consider extended logical space even though the part of the y+yShift curve
   -- that is not in this cell does not get used in the integral because the fixed
   -- y limit cuts if off. But if we didn't consider extended logical space the shape of the
   -- x-limit function would not come out right.
--   print("ySN= ",ySNIn, " | y0= ",y0In)
   local eps = tol  -- Need wiggle room because it might fail if the root is at the boundary.
   local xOutN = root.ridders(lossF, xLim[1]-eps, xLim[2]+eps, rootStop(tol))
   return xOutN
end


-- Find if a table of tables has a nil.
local nilAny2x2 = function(tblIn)
  for i = 1,2 do
     for j = 1,2 do
        if tblIn[i][j]==nil then return true end
     end
  end
  return false
end
-- Count how many nils there are in a table of tables.
local nilHowmany2x2 = function(tblIn)
  local count, nilIdxs, nonNilIdxs = 0, {}, {}
  for i = 1,2 do
     for j = 1,2 do
        if tblIn[i][j] == nil then
           count = count+1
           table.insert(nilIdxs,{i,j})
        else
           table.insert(nonNilIdxs,{i,j})
        end
     end
  end
  return count, nilIdxs, nonNilIdxs
end

-- Determine the cells that each cell in fldDest needs.
local numBasis1D = basis1D:numBasis()
local basis1Dev  = Lin.Vec(numBasis1D)

local idx1D      = Lin.IntVec(1)
local indexer    = fldSrc:genIndexer()
local indexer1D  = yShift:genIndexer()
local localRange = fldSrc:localRange()
local srcPtr, destPtr, yShiftPtr = fldSrc:get(1), fldDest:get(1), yShift:get(1)
local xc, dx     = Lin.Vec(2), Lin.Vec(2)
local xcS, dxS   = Lin.Vec(2), Lin.Vec(2)
local stepFrac   = {0.1,0.1}          -- Size of step taken around the boundary, as a fraction of cell length.
local srcCells   = {}
local numSteps   = {math.ceil(1/stepFrac[1]),math.ceil(1/stepFrac[2])}

-- Set up DG projection like EvalOnNodes:
local nodes    = SerendipityNodes["nodes1xp" .. polyOrder1D]
local numNodes = #nodes
local numVal   = 1
--local fvLo, fvUp = Lin.Mat(numNodes, numVal), Lin.Mat(numNodes, numVal)
local fv = Lin.Mat(numNodes, numVal)
local compToPhysTempl = xsys.template([[
return function (eta, dx, xc, xOut)
|for i = 1, NDIM do
   xOut[${i}] = 0.5*dx[${i}]*eta[${i}] + xc[${i}]
|end
end
]])
local evalFuncTempl = xsys.template([[
return function (tCurr, xn, func, fout)
|for i = 1, M-1 do
   fout[${i}],
|end
   fout[${M}] = func(tCurr, xn)
end
]])
ffi.cdef[[
  void nodToMod(double* fN, int numNodes, int numVal, int ndim, int p, double* fM);
]]
local evalFunc = loadstring(evalFuncTempl { M = numVal } )()
local qL_y, qU_y = Lin.Vec(numNodes), Lin.Vec(numNodes)
local compToPhys = loadstring(compToPhysTempl {NDIM = 1} )()
local projXlim = function(funcIn, q_yOut)
   for i = 1, numNodes do
      local xj = nodes[i]
      evalFunc(0.0, xj, funcIn, fv[i])
   end
   ffiC.nodToMod(fv:data(), numNodes, numVal, 1, polyOrder1D, q_yOut:data())
end
local xi = Lin.Vec(1)                  -- Coordinates at node.
local dyPartial, ycPartial = Lin.Vec(1), Lin.Vec(1)
local projXlimGen = function(funcIn, dxIn, xcIn, q_yOut)
   for i = 1, numNodes do
      compToPhys(nodes[i], dxIn, xcIn, xi)      -- Compute physical coordinate xi.
      evalFunc(0.0, xi, funcIn, fv[i])
   end
   ffiC.nodToMod(fv:data(), numNodes, numVal, 1, polyOrder1D, q_yOut:data())
end

local interp = TwistShiftInterpDecl.selectTwistShiftInterp(basis:id(), 2, polyOrder1D)
local interpGen = TwistShiftInterpDecl.selectTwistShiftInterpGen(basis:id(), 2, polyOrder1D)
local interpGenSub = TwistShiftInterpDecl.selectTwistShiftInterpGenSub(basis:id(), 2, polyOrder1D)

--if (1==2) then

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
   fldDest:fill(indexer(idx), destPtr)

   local cellLower  = {grid:cellLowerInDir(1),grid:cellLowerInDir(2)}
   local cellUpper  = {grid:cellUpperInDir(1),grid:cellUpperInDir(2)}
   local cellBounds = { {cellLower[1],cellUpper[1]}, {cellLower[2],cellUpper[2]} }
   local stepSize   = {grid:dx(1)/numSteps[1],grid:dx(2)/numSteps[2]}

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

   print(string.format("idx = (%d,%d)",idx[1],idx[2]))

--   if idx[1]==1 and idx[2]==1 then
                                  
   for iC = 1, #srcCells do
      print(string.format("   from = (%d,%d)",srcCells[iC][1],srcCells[iC][2]))
      -- In each contributing cell, approximate the functions qInvL(y) and qInvU(y)
      -- for the min/max of the x-integral with a polynomial. Compute the coefficients
      -- of that polynomial with a projection of yShiftNormInvFunc onto the local basis.
      local idxS = srcCells[iC]

      grid:setIndex(idxS)
      grid:cellCenter(xcS)
      grid:getDx(dxS)
      fldSrc:fill(indexer(idxS), srcPtr)
      local srcCellBounds = { {grid:cellLowerInDir(1),grid:cellUpperInDir(1)},
                              {grid:cellLowerInDir(2),grid:cellUpperInDir(2)} }

--      if (idxS[1]==1 and idxS[2]==6)       -- the following is temporary code to test the case of a contribution in which
                                           -- the y+yShift curves do not cross the x-boundaries of the cell.
--        or (idxS[1]==1 and idxS[2]==10)       -- the following is temporary code to test the case of a contribution in which
--      then

      -- Find the four points where yLowerDest+yShift yUpperDest+yShift intersect
      -- the yLowerSrc and yUpperSrc lines.
      local trapCornerX = {{nil, nil}, {nil, nil}}
      for i = 1, 2 do
         for j = 1, 2 do 
            local yDest = cellBounds[2][i] 
            local ySrc  = srcCellBounds[2][j] 
            trapCornerX[i][j] = yShiftedyIntersect(cellBounds[1], yDest, ySrc)
--            if trapCornerX[i][j]==nil then
--               print(string.format("trapCornerX[%d][%d] = ",i,j),nil)
--            else
--               print(string.format("trapCornerX[%d][%d] = %f",i,j,trapCornerX[i][j]))
--            end
         end
      end

      if not nilAny2x2(trapCornerX) then
         -- Scenario sN.
         -- The contributing region in the source cell is 4-sided with all vertices
         -- contained within the source cell.

         -- Establish if the yLowerDest+yShift/yUpperDest+yShift defines the lower
         -- or upper limit of the x-integral below. Do so by comparing where y=yLowerSrc
         -- crosses yLowerDest+yShift and yUpperDest+yShift.
         local y0Lo, y0Up, yLimLo, yLimUp
         local xLimLoLo, xLimLoUp, xLimUpLo, xLimUpUp
         if trapCornerX[1][1] < trapCornerX[2][1] then
            -- Monotonically decreasing yShift.
            y0Lo, y0Up = cellBounds[2][1], cellBounds[2][2]
            xLimLoLo, xLimLoUp = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
            xLimUpLo, xLimUpUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
         else
            -- Monotonically increasing yShift.
            y0Lo, y0Up = cellBounds[2][2], cellBounds[2][1]
            xLimLoLo, xLimLoUp = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
            xLimUpLo, xLimUpUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
         end

         -- Project the function describing the limits of the x-integral onto the basis.
         -- This projection is modeled after Update/EvalOnNodes.lua.
--         local evaluateLo = function(tCurr,xn) return yShiftNormInvFunc(xn[1], y0Lo, xcS, xc, dx) end
--         local evaluateUp = function(tCurr,xn) return yShiftNormInvFunc(xn[1], y0Up, xcS, xc, dx) end
         local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
         local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
         --projXlim(evaluateLo, qL_y)
         --projXlim(evaluateUp, qU_y)
         yLimLo, yLimUp = -1., 1.
         dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
         projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
         projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
         --for i = 1, numNodes do
         --   print(qL_y[i], qU_y[i])
         --end
--         for k=1,numBasis1D do
--            print(string.format("qL[%d], qU[%d] = %f, %f",k-1, k-1, qL_y:data()[k-1], qU_y:data()[k-1]))
--         end

         -- Having obtained a polynomial expansion to the limits of the x-integral,
         -- proceed to add the contribution from this source cell.
         --interp(qL_y:data(), qU_y:data(), srcPtr:data(), destPtr:data())
         interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())

         --for i=1,4 do
         --   print(string.format("src[%d] = %f | dest[%d] = %f",i-1,srcPtr:data()[i-1],i-1,destPtr:data()[i-1]))
         --end
      --end

      --elseif (idxS[1]==1 and idxS[2]==8)   -- the following is temporary code to test the case of a contribution in which
      --                                     -- only one y+yShift curve intersects the source cell.
      --      or (idxS[1]==1 and idxS[2]==7)   -- the following is temporary code to test the case of a contribution in which
      --                                       -- one y+yShift curve intersects one y-boundary of the source cell and
      --                                       -- the other y+yShift curve intersects both y-boundaries.
      --      then

      ---- Find the four points where yLowerDest+yShift yUpperDest+yShift intersect
      ---- the yLowerSrc and yUpperSrc lines.
      --local trapCornerX = {{nil, nil}, {nil, nil}}
      --for i = 1, 2 do   -- Loop over yDest lower/upper.
      --   for j = 1, 2 do    -- Loop over ySrc lower/upper.
      --      local yDest = cellBounds[2][i] 
      --      local ySrc  = srcCellBounds[2][j] 
      --      trapCornerX[i][j] = yShiftedyIntersect(cellBounds[1], yDest, ySrc)
      --   end
      --end

      else

      local nilNum, nilIdxs, nonNilIdxs = nilHowmany2x2(trapCornerX)
      
      if nilNum == 3 then
         -- There are 4 possible scenarios:
         --   si)   yLowerDest+yShift intersects the lower x-boundary.
         --   sii)  yLowerDest+yShift intersects the upper x-boundary.
         --   siii) yUpperDest+yShift intersects the lower x-boundary.
         --   siv)  yUpperDest+yShift intersects the upper x-boundary.
         -- Define yL/yU to be the lower/upper limit of the y-integral.
         -- Scenarios a,b have an yU=1 and TBD yL. Scenarios c,d have yL=-1 and TBD yU.
         -- Scenarios a,c have xL=-1 and variable xU. Scenarios b,d have xU=1 and variably xL.

         -- Establish if the yLowerDest+yShift/yUpperDest+yShift defines the lower
         -- or upper limit of the x-integral below. Do so by comparing where y=yLowerSrc
         -- crosses yLowerDest+yShift and yUpperDest+yShift.
         local y0Lo, y0Up, yLimLo, yLimUp
         local evaluateLo, evaluateUp
         if nonNilIdxs[1][1]==1 then
            -- Scenario si or sii.
            yLimUp = 1.
            local yShiftUp = yShiftFunc(srcCellBounds[1][2])
            y0Lo = cellBounds[2][1]
            if yShiftUp >= yShiftFunc(trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]]) then
               yLimLo = wrapNum(cellLower[2] + yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])  -- Scenario si.
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               xLimLo = -1
               xLimUp = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return -1.0 end
               --evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
            else
               yLimLo = wrapNum(cellLower[2] + yShiftUp,grid:lower(2),grid:upper(2),idxS[2]==numCells[2])   -- Scenario sii.
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               xLimLo = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               xLimUp = 1
               --evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
            end
         else
            -- Scenario siii or siv.
            yLimLo = -1.
            local yShiftUp = yShiftFunc(srcCellBounds[1][2])
            y0Up = cellBounds[2][2]
            if yShiftUp <= yShiftFunc(trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]]) then
               yLimUp = wrapNum(cellUpper[2]+yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])   -- Scenario siii.
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               xLimLo = -1
               xLimUp = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return -1.0 end
               --evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
            else
               yLimUp = wrapNum(cellUpper[2] + yShiftUp,grid:lower(2),grid:upper(2),idxS[2]==numCells[2])   -- Scenario siv.
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               xLimLo = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               xLimUp = 1
               --evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
            end
         end
         dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
         projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
         --projXlim(evaluateUp, qU_y)
         projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
         --print(string.format("xLimLoLo=%f | xLimLoUp=%f | xLimUpLo=%f | xLimUpUp=%f | yLimLo=%f | yLimUp=%f",xLimLoLo,xLimLoUp,xLimUpLo,xLimUpUp,yLimLo,yLimUp))
--         print(string.format("xLimLo=%f | xLimUp=%f | yLimLo=%f | yLimUp=%f",xLimLo,xLimUp,yLimLo,yLimUp))
--         for i=1,numBasis1D do
--            print(string.format("qL[%d]=%f | qU[d]=%f",i-1,qL_y:data()[i-1],i-1,qU_y:data()[i-1]))
--         end


         -- Having obtained a polynomial expansion to the limits of the x-integral and the
         -- values of the y-integral limits, add the contribution from this source cell.
         interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())

      elseif nilNum == 1 then
         -- There are 4 possible scenarios:
         --   sv)    yUpperDest+yShift doesn't intersect the lower y-boundary and intersects the lower x-boundary.
         --   svi)   yUpperDest+yShift doesn't intersect the lower y-boundary and intersects the upper x-boundary..
         --   svii)  yLowerDest+yShift doesn't intersect the upper y-boundary and intersects the lower x-boundary.
         --   sviii) yLowerDest+yShift doesn't intersect the upper y-boundary and intersects the upper x-boundary.

         -- Split the contribution to this cell into a contribution like those in scenarios si-siv above, and
         -- a contribution like scenario sN (4-sided region whole contained).

         if (nilIdxs[1][1]==2) and (nilIdxs[1][2]==1) then
            -- Scenario sv or svi.

            -- Establish if the yLowerDest+yShift/yUpperDest+yShift defines the lower
            -- or upper limit of the x-integral below. Do so by comparing where y=yUpperSrc
            -- crosses yLowerDest+yShift and yUpperDest+yShift.
            local y0Lo, y0Up, yLimLo, yLimUp
            local xLimLoLo, xLimLoUp, xLimUpLo, xLimUpUp
            if trapCornerX[1][2] > trapCornerX[2][2] then
               -- Scenario sv.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][2], cellBounds[2][1]
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               projXlim(evaluateLo, qL_y)
               projXlim(evaluateUp, qU_y)
               yLimLo = wrapNum(cellUpper[2] + yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               yLimUp = 1.
               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               -- Add scenario siii-like contribution.
               y0Up = cellBounds[2][1]
               yLimUp = yLimLo    -- Keep the order of yLimUp and yLimLo here.
               yLimLo = -1.
               evaluateLo = function(tCurr,xn) return -1.0 end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               projXlim(evaluateLo, qL_y)
               projXlim(evaluateUp, qU_y)
               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
            else
               -- Scenario svi.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][1], cellBounds[2][2]
--               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
--               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
--               projXlim(evaluateLo, qL_y)
--               projXlim(evaluateUp, qU_y)
               xLimLoLo, xLimLoUp = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               xLimUpLo, xLimUpUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), 1.
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               yLimLo = wrapNum(cellUpper[2] + yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               yLimUp = 1.
--               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
               interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario siv-like contribution.
               yLimUp = yLimLo    -- Keep the order of yLimUp and yLimLo here.
               yLimLo = -1.
               y0Lo   = cellBounds[2][1]
--               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               xLimLo = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
--               projXlim(evaluateLo, qL_y)
--               projXlim(evaluateUp, qU_y)
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
--               print("here")
--               print(string.format("xLimLo=%f | xLimUp=%f | yLimLo=%f | yLimUp=%f",xLimLo,xLimUp,yLimLo,yLimUp))
--               print(string.format("dyP=%f | ycP=%f",dyPartial[1],ycPartial[1]))
               projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
--               print("here too")
               projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
--               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
            end

         else
            -- Scenario svii or sviii.
--            for i=1,2 do
--               for j=1,2 do
--                  if trapCornerX[i][j]==nil then
--                     print(string.format("trapCornerX[%d][%d] = ",i,j),nil)
--                  else
--                     print(string.format("trapCornerX[%d][%d] = %f",i,j,trapCornerX[i][j]))
--                  end
--               end
--            end

            -- Establish if the yLowerDest+yShift/yUpperDest+yShift defines the lower
            -- or upper limit of the x-integral below. Do so by comparing where y=yLowerSrc
            -- crosses yLowerDest+yShift and yUpperDest+yShift.
            local y0Lo, y0Up, yLimLo, yLimUp
            local xLimLoLo, xLimLoUp, xLimUpLo, xLimUpUp
            if trapCornerX[1][1] < trapCornerX[2][1] then
               -- Scenario svii.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][1], cellBounds[2][2]
               xLimLoLo, xLimLoUp = -1., (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               xLimUpLo, xLimUpUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
--               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
--               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               --projXlim(evaluateLo, qL_y)
               --projXlim(evaluateUp, qU_y)
               yLimLo = -1.
               yLimUp = wrapNum(cellLower[2] + yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
--               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
               interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario si-like contribution.
               yLimLo = yLimUp    -- Keep the order of yLimUp and yLimLo here.
               yLimUp = 1.
               xLimLo = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
               y0Up = cellBounds[2][2]
               evaluateLo = function(tCurr,xn) return -1.0 end
--               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
--               projXlim(evaluateLo, qL_y)
--               projXlim(evaluateUp, qU_y)
--               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qL_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qU_y)
               interpGenSub(qL_y:data(), qU_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
            else
               -- Scenario sviii.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][2], cellBounds[2][1]
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Up, xcS, xc, dx) end
               projXlim(evaluateLo, qL_y)
               projXlim(evaluateUp, qU_y)
               yLimLo = -1.
               yLimUp = wrapNum(cellLower[2] + yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
               -- Add scenario sii-like contribution.
               yLimLo = yLimUp    -- Keep the order of yLimUp and yLimLo here.
               yLimUp = 1.
               y0Lo   = cellBounds[2][2]
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncEx(xn[1], y0Lo, xcS, xc, dx) end
               evaluateUp = function(tCurr,xn) return 1.0 end
               projXlim(evaluateLo, qL_y)
               projXlim(evaluateUp, qU_y)
               interpGen(qL_y:data(), qU_y:data(), yLimLo, yLimUp, srcPtr:data(), destPtr:data())
            end

         end

      end
--      end

   end
--   end

end
end

fldDest:write("fldDest.bp")

local intFldDest = DataStruct.DynVector {
   numComponents = 1,
}
intQuant:advance(0., {fldDest}, {intFldDest})
intFldDest:write("intFldDest.bp", 0., 0)
