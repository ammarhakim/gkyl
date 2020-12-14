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
local math = require "sci.math"  -- For sign function.

-- Create two fields on a 2D grid. Then interpolate one field onto the other
-- but with a shift in y that is a function of x (assuming periodicity in y).

local polyOrder    = 1
local lower        = {-2.0, -1.50}
local upper        = { 2.0,  1.50}
local numCells     = {36, 30}
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
-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(x)
--                      return 1./(1.+0.25*x)
                      return 0.3*x+2.
--                      return 0.1*(x+2.)^2+0.2*x+0.6
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
                 local muX, muY   = 0., 0.
                 local sigX, sigY = 0.3, 0.3
                 return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
--                 return 1.
--                 return yShiftFunc(x)*math.sin((2.*math.pi/(upper[2]-lower[2]))*y)
--                  return y-yShiftFunc(x) --+2  --yShiftInvFunc(yShiftFunc(x))
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

-- Find the root of the yShiftFunc = 0 equation. If it is found it means
-- that yShiftFunc needs to be changed to something that doesn't go through zero.
local function checkForRoot(funcIn, xBounds, tol) 
   local lossFl, lossFu = funcIn(xBounds[1]), funcIn(xBounds[2])
   if math.abs(lossFl) < tol  then
      return xBounds[1]
   elseif math.abs(lossFu) < tol  then
      return xBounds[2]
   else
      if lossFl*lossFu < 0. then
         return root.ridders(funcIn, xBounds[1], xBounds[2], rootStop(tol))
      else
         return nil
      end
   end
end
xRoot = checkForRoot(yShiftFunc, {grid:lower(1), grid:upper(1)}, 1.e-13)
if xRoot then
   assert( xRoot==nil, string.format("yShiftFunc is zero at %f. Please select a function that is not zero in the domain.",xRoot))
end


local function sign(num)
--   if math.abs(num) < 1.e-14 then
--      return 0
--   elseif num > 0 then
   if num > 0 then
      return 1
   elseif num < 0 then
      return -1
   else
      return 0
   end
end

-- Given a y-coordinate of the destination cell (yDest), and a y-coordinate of the
-- source cell (ySrc), find the x-coordinate of the point where the y+yShift(x) and
-- y=ySrc lines intersect. Need to provide a window in x where to search ([xBounds[1],xBounds[2]]).
local yShiftedyIntersect = function(xBounds, yDest, ySrc, jDest, jSrc)
   local Ly = grid:upper(2)-grid:lower(2)
   local function lossF(xIn)
      -- Rather than wrapping yDest+yShift, it might be more robust
      -- to look for the intersection in an extended space.
      local shift     = yShiftFunc(xIn)
      local shiftSign = sign(shift)
      if math.abs(shift) > (Ly+1.e-10) then
      print("xIn=", xIn, "| ySrc=",ySrc, " | yDest=",yDest," | shift=",shift)
         shift = shiftSign*wrapNum(math.abs(shift), 0., Ly, true)
      end
      local shiftedY = yDest+shift
      if shiftSign==1 and ySrc < yDest then
         return shiftedY - (grid:upper(2) + ySrc-grid:lower(2))
      elseif shiftSign==-1 and ySrc > yDest then
         return shiftedY - (grid:lower(2) - (grid:upper(2)-ySrc))
      else
      print("xIn=", xIn, "| ySrc=",ySrc, " | yDest=",yDest," | shift=",shift," | diff=",shiftedY - ySrc)
         return shiftedY - ySrc 
      end
--      if shiftedY > grid:upper(2) then
----         print("x-xLower, x-xUpper =", xIn-xBounds[1], xIn-xBounds[2])
--         print("xIn=", xIn," | shiftedY=",shiftedY-1.5, "| ySrc=",ySrc, " | yDest=",yDest, " loss=",(grid:upper(2) + ySrc-grid:lower(2)))
--         return shiftedY - (grid:upper(2) + ySrc-grid:lower(2))
----         return shiftedY - (grid:upper(2) + ySrc-grid:lower(2))
----         return wrapNum(shiftedY,grid:lower(2),grid:upper(2),true) - ySrc
--      elseif shiftedY < grid:lower(2) then
--         return shiftedY - (grid:lower(2) - (grid:upper(2)-ySrc))
--      else
--         return shiftedY - ySrc 
--      end
--      if math.abs(shiftedY-ySrc) > 2.1*grid:dx(2) then
--         print("abs(shiftedY-ySrc)=",math.abs(shiftedY-ySrc)," | dy=",grid:dx(2))
--         print("shiftedY=",shiftedY," | wrap=",wrapNum(shiftedY,grid:lower(2),grid:upper(2),true)," | ySrc=",ySrc, " | loss=",wrapNum(shiftedY,grid:lower(2),grid:upper(2),true) - ySrc)
--         return wrapNum(shiftedY,grid:lower(2),grid:upper(2),true) - ySrc
--      else
--         return shiftedY - ySrc 
--      end
--      return shiftedY - ySrc 
   end
   local function lossF2(xIn)
      local shiftedY = yDest+yShiftFunc(xIn)
      return wrapNum(shiftedY,grid:lower(2),grid:upper(2),true) - ySrc
   end
   local tol = 1.e-13
   local eps = 0. --1.e-14
   -- Check if the root is bracketed by xBounds (otherwise root.ridders throws an error).
   -- We however use the lack of a root to identify the topology of the source cell region.
   local function rootFind(lossFunc)
      local lossFl, lossFu = lossFunc(xBounds[1]+eps), lossFunc(xBounds[2]-eps)
   print("xBounds[1], xBounds[2] = ",xBounds[1], xBounds[2])
   print("lossFl, lossFu = ",lossFl, lossFu)
--      if lossFl*lossFu < 0. then
--         return root.ridders(lossFunc, xBounds[1]-eps, xBounds[2]+eps, rootStop(tol))
--      else
--         if math.abs(lossFl) < tol  then
--            return xBounds[1]
--         elseif math.abs(lossFu) < tol  then
--            return xBounds[2]
--         else
--            return nil
--         end
--      end
      if math.abs(lossFl) < tol  then
         return xBounds[1]
      elseif math.abs(lossFu) < tol  then
         return xBounds[2]
      else
         if lossFl*lossFu < 0. then
            return root.ridders(lossFunc, xBounds[1]-eps, xBounds[2]+eps, rootStop(tol))
         else
            -- It's possible that the y+yShift line crosses the domain y-boundary.
            -- Need a different loss function that accounts for that.
            local yShLo, yShUp = yDest+yShiftFunc(xBounds[1]), yDest+yShiftFunc(xBounds[2])
            return nil
         end
      end
   end
   return rootFind(lossF)
   --local root0 = rootFind(lossF)
   --if root0 == nil then
   --   -- It's possible that the y+yShift intersects a periodic copy of this domain.
   --   -- Wrap the y+yShift curve back to this domain and try again.
   --   return rootFind(lossF2)
   --else
   --   return root0
   --end
end

-- Given a y-coordinate of the destination cell (yDest), and a y-coordinate of the
-- source cell (ySrc), find the x-coordinate of the point where the y+yShift(x) and
-- y=ySrc lines intersect. Need to provide a window in x where to search ([xBounds[1],xBounds[2]]).
-- If y+yShift-ySrc=0 has no roots, it is possible that y+yShift intersects a periodic
-- copy of this domain. Check for such cases by looking for the roots of
-- y+yShift-(ySrc+N*Ly)=0 where Ly is the length of the domain along y and N is an integer.
local yShiftedyIntersectEx = function(xBounds, yDest, ySrc, jDest, jSrc)
   local Ly = grid:upper(2)-grid:lower(2)
   local function lossF(xIn)
      return yDest+yShiftFunc(xIn) - ySrc 
   end
   local tol = 1.e-13
   local function rootFind(lossFunc)
       -- Check if the root is bracketed by xBounds (otherwise root.ridders throws an error).
      local lossFl, lossFu = lossFunc(xBounds[1]), lossFunc(xBounds[2])
      if math.abs(lossFl) < tol  then
         return xBounds[1]
      elseif math.abs(lossFu) < tol  then
         return xBounds[2]
      else
         if lossFl*lossFu < 0. then
            return root.ridders(lossFunc, xBounds[1], xBounds[2], rootStop(tol))
         else
            return nil
         end
      end
   end
   root0 = rootFind(lossF)
   if root0 == nil then
      -- Maybe the y+yShift intersects y=ySrc in a periodic copy of this domain.
      -- Find the roots of y+yShift-(ySrc+nP*Ly)=0 where Ly is the length of the
      -- domain along y and nP is an integer.
      local nP = {} 
      -- Evaluate the y+yShift at some points in [xBounds[1],xBounds[2]] and obtain
      -- potential nP multipliers.
      local stepFrac = 0.1
      local numSteps = math.ceil(1./stepFrac)
      local stepSize = (xBounds[2]-xBounds[1])/numSteps
      for sI = 0, numSteps do
         local xp      = xBounds[1]+sI*stepSize 
         local exShift = (yDest+yShiftFunc(xp)-grid:lower(2))/Ly
         local newN    = math.sign(exShift)*math.floor(math.abs(exShift))
--         print(string.format("xp=%f | yDest+yShift=%f | newN=%d",xp,(yDest+yShiftFunc(xp)),newN))
         if lume.find(nP,newN)==nil then table.insert(nP,newN) end
      end
      for _, iN in ipairs(nP) do
        local function lossFex(xIn)
           return yDest+yShiftFunc(xIn) - (ySrc+iN*Ly)
        end
        root0 = rootFind(lossFex)
        if root0 then break end
      end
   end
   return root0
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
   local tol = 1.e-11
   -- Have to consider extended logical space even though the part of the y+yShift curve
   -- that is not in this cell does not get used in the integral because the fixed
   -- y limit cuts if off. But if we didn't consider extended logical space the shape of the
   -- x-limit function would not come out right.
--   print("ySN= ",ySNIn, " | y0= ",y0In)
   local eps = tol  -- Need wiggle room because it might fail if the root is at the boundary.
   local xOutN = root.ridders(lossF, xLim[1]-eps, xLim[2]+eps, rootStop(tol))
   return xOutN
end


-- Count how many nils there are in a table.
local nilHowmany2 = function(tblIn)
  local count, nilIdxs, nonNilIdxs = 0, {}, {}
  for i = 1,2 do
     if tblIn[i] == nil then
        count = count+1
        table.insert(nilIdxs,i)
     else
        table.insert(nonNilIdxs,i)
     end
  end
  return count, nilIdxs, nonNilIdxs
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
local qLo_x, qUp_x = Lin.Vec(numNodes), Lin.Vec(numNodes)
local qLo_y, qUp_y = Lin.Vec(numNodes), Lin.Vec(numNodes)
local qLoP_y, qUpP_y = Lin.Vec(numNodes), Lin.Vec(numNodes)
local compToPhys = loadstring(compToPhysTempl {NDIM = 1} )()
local projXlim = function(funcIn, q_yOut)
   for i = 1, numNodes do
      local xj = nodes[i]
      evalFunc(0.0, xj, funcIn, fv[i])
   end
   ffiC.nodToMod(fv:data(), numNodes, numVal, 1, polyOrder1D, q_yOut:data())
end
local xi = Lin.Vec(1)                  -- Coordinates at node.
local dxPartial, xcPartial = Lin.Vec(1), Lin.Vec(1)
local dyPartial, ycPartial = Lin.Vec(1), Lin.Vec(1)
local dyPartialP, ycPartialP = Lin.Vec(1), Lin.Vec(1)
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
local interpGenSubY = TwistShiftInterpDecl.selectTwistShiftInterpGenSubY(basis:id(), 2, polyOrder1D)
local interpGenSubMtwoCorners = TwistShiftInterpDecl.selectTwistShiftInterpGenSubMinusTwoCorners(basis:id(), 2, polyOrder1D)

--if (1==2) then

for idx in localRange:rowMajorIter() do

   -- One way of determining which cells this cell is magnetically connected to
   -- is to take points around the boundary of the cell and calculate the shift from those
   -- points, then find out which cells own the shifted points. This is perhaps limited to
   -- grid that are not refined in the z-direction; in that case we would have to consider
   -- more than just cell-boundary points.
   -- For now simply take small steps around the boundaries of the 2D cell, compute
   -- the shifted points and their owner at each step.

--   if idx[1]==36 and idx[2]==6 then

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
   local delta    = {1.e-9*grid:dx(1),1.e-9*grid:dx(2)}
   for dC = 1,2 do            -- Loop over x=const and y=const boundaries.
      for xS = 1,2 do         -- Loop over lower/upper boundary.
--         print("dC=",dC," | xS=",xS)

         -- Do the first point separately (it requires pickLower=false in findCell).
         local evPoint = {cellLower[1]+delta[1], cellLower[2]+delta[2]}
         evPoint[dC]   = evPoint[dC]+(xS-1)*(grid:dx(dC)-2.*delta[dC])

         basis1D:evalBasis({(evPoint[1]-grid:cellCenterInDir(1))/(0.5*grid:dx(1))}, basis1Dev)
         local yShiftB = 0.   -- y-shift at boundary.
         for k = 1, numBasis1D do
            yShiftB = yShiftB + yShiftPtr[k]*basis1Dev[k]
         end

--         print(string.format("yShifted = %f | idx[1]=%d",evPoint[2]+yShiftB,idx[1]))
--         print(string.format("Find point=(%f,%f)",evPoint[1],wrapNum(evPoint[2]+yShiftB,grid:lower(2),grid:upper(2),true)))
         local idxShifted0   = {nil,nil}
         local chooseLower   = (dC==2 and xS==2) and true or false
         local wrapInclusive = (dC==2 and xS==2) and true or false
         local searchPoint   = {evPoint[1],wrapNum(evPoint[2]+yShiftB,grid:lower(2),grid:upper(2),true)}
         grid:findCell(searchPoint,idxShifted0,chooseLower,{idx[1],nil})
         if lume.findTable(srcCells,idxShifted0)==nil then table.insert(srcCells,idxShifted0) end
         
         local stepDim, newP = dC % 2+1, {nil,nil}
         newP[dC] = evPoint[dC]
         for sI = 1, numSteps[stepDim] do

            newP[stepDim] = evPoint[stepDim] + sI*stepSize[stepDim]
            if sI==numSteps[stepDim] then newP[stepDim]=newP[stepDim]-2.*delta[stepDim] end
         
            -- Evaluate yShift at this point.
            basis1D:evalBasis({(newP[1]-grid:cellCenterInDir(1))/(0.5*grid:dx(1))}, basis1Dev)
            local yShiftB = 0.   -- y-shift at boundary.
            for k = 1, numBasis1D do
               yShiftB = yShiftB + yShiftPtr[k]*basis1Dev[k]
            end

            -- Find the index of the cell that owns the shifted point.
            local idxShifted  = {nil,nil}
--         print(string.format("Find point=(%f,%f)",newP[1],wrapNum(newP[2]+yShiftB,grid:lower(2),grid:upper(2),true)))
            local searchPoint = {newP[1],wrapNum(newP[2]+yShiftB,grid:lower(2),grid:upper(2),true)}
            grid:findCell(searchPoint,idxShifted,true,{idx[1],nil})
            if lume.findTable(srcCells,idxShifted)==nil then table.insert(srcCells,idxShifted) end
         end
      end
   end

   print(string.format("idx = (%d,%d)",idx[1],idx[2]))

--   if idx[1]==1 and idx[2]==1 then
                                  
   for iC = 1, #srcCells do
      -- In each contributing cell, approximate the functions qInvL(y) and qInvU(y)
      -- for the min/max of the x-integral with a polynomial. Compute the coefficients
      -- of that polynomial with a projection of yShiftNormInvFunc onto the local basis.
      local idxS = srcCells[iC]
      print(string.format("   from = (%d,%d)",idxS[1],idxS[2]))

                                  
      grid:setIndex(idxS)
      grid:cellCenter(xcS)
      grid:getDx(dxS)
      fldSrc:fill(indexer(idxS), srcPtr)
      local srcCellBounds = { {grid:cellLowerInDir(1),grid:cellUpperInDir(1)},
                              {grid:cellLowerInDir(2),grid:cellUpperInDir(2)} }

--      if (idxS[1]==3 and idxS[2]<4)       -- the following is temporary code to test the case of a contribution in which
                                           -- the y+yShift curves do not cross the x-boundaries of the cell.
--        or (idxS[1]==1 and idxS[2]==10)       -- the following is temporary code to test the case of a contribution in which
--      then

      -- Find the four points where yLowerDest+yShift and yUpperDest+yShift
      -- intersect the yLowerSrc and yUpperSrc lines.
      local trapCornerX = {{nil, nil}, {nil, nil}}
      for i = 1, 2 do
         for j = 1, 2 do 
            local yDest = cellBounds[2][i] 
            local ySrc  = srcCellBounds[2][j] 
            trapCornerX[i][j] = yShiftedyIntersectEx(cellBounds[1], yDest, ySrc, idx[2], idxS[2])
            if trapCornerX[i][j]==nil then
               print(string.format("trapCornerX[%d][%d] = ",i,j),nil)
            else
               print(string.format("trapCornerX[%d][%d] = %f",i,j,trapCornerX[i][j]))
            end
         end
      end

--   if idx[1]==1 and idx[2]<-9 then
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
         local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
         local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
         yLimLo, yLimUp = -1., 1.
         dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
         projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
         projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)

         -- Having obtained a polynomial expansion to the limits of the x-integral,
         -- proceed to add the contribution from this source cell.
         interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())

      --end

      --elseif (idxS[1]==1 and idxS[2]==8)   -- the following is temporary code to test the case of a contribution in which
      --                                     -- only one y+yShift curve intersects the source cell.
      --      or (idxS[1]==1 and idxS[2]==7)   -- the following is temporary code to test the case of a contribution in which
      --                                       -- one y+yShift curve intersects one y-boundary of the source cell and
      --                                       -- the other y+yShift curve intersects both y-boundaries.
      --      then

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
               local xLimLo = -1
               local xLimUp = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return -1.0 end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
            else
               yLimLo = wrapNum(cellLower[2] + yShiftUp,grid:lower(2),grid:upper(2),idxS[2]==numCells[2])   -- Scenario sii.
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               local xLimLo = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               local xLimUp = 1
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
               local xLimLo = -1
               local xLimUp = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return -1.0 end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
            else
               yLimUp = wrapNum(cellUpper[2] + yShiftUp,grid:lower(2),grid:upper(2),idxS[2]==numCells[2])   -- Scenario siv.
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               local xLimLo = (trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]] - xcS[1])/(0.5*dxS[1])
               local xLimUp = 1
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
            end
         end
         dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
         projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
         projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
--         --print(string.format("xLimLoLo=%f | xLimLoUp=%f | xLimUpLo=%f | xLimUpUp=%f | yLimLo=%f | yLimUp=%f",xLimLoLo,xLimLoUp,xLimUpLo,xLimUpUp,yLimLo,yLimUp))
--         print(string.format("xLimLo=%f | xLimUp=%f | yLimLo=%f | yLimUp=%f",xLimLo,xLimUp,yLimLo,yLimUp))
--         for i=1,numBasis1D do
--            print(string.format("qL[%d]=%f | qU[d]=%f",i-1,qLo_y:data()[i-1],i-1,qUp_y:data()[i-1]))
--         end


         -- Having obtained a polynomial expansion to the limits of the x-integral and the
         -- values of the y-integral limits, add the contribution from this source cell.
         interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
--         for i=1,basis:numBasis() do
--            print(string.format("src[%d]=%f | dest[%d]=%f",i-1,srcPtr:data()[i-1],i-1,destPtr:data()[i-1]))
--         end

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
               xLimLoLo, xLimLoUp = -1., (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               xLimUpLo, xLimUpUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               yLimLo = wrapNum(cellUpper[2] + yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               yLimUp = 1.
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario siii-like contribution.
               yLimUp = yLimLo    -- Keep the order of yLimUp and yLimLo here.
               yLimLo = -1.
               y0Up = cellBounds[2][1]
               xLimLo = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return -1.0 end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               print("here")
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
            else
               -- Scenario svi.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][1], cellBounds[2][2]
               xLimLoLo, xLimLoUp = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               xLimUpLo, xLimUpUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), 1.
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               yLimLo = wrapNum(cellUpper[2] + yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               yLimUp = 1.
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario siv-like contribution.
               yLimUp = yLimLo    -- Keep the order of yLimUp and yLimLo here.
               yLimLo = -1.
               y0Lo   = cellBounds[2][1]
               xLimLo = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
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
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               yLimLo = -1.
               yLimUp = wrapNum(cellLower[2] + yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario si-like contribution.
               yLimLo = yLimUp    -- Keep the order of yLimUp and yLimLo here.
               yLimUp = 1.
               xLimLo = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
               y0Up = cellBounds[2][2]
               evaluateLo = function(tCurr,xn) return -1.0 end
               evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
            else
               -- Scenario sviii.
               -- Add scenario sN-like contribution.
               y0Lo, y0Up = cellBounds[2][2], cellBounds[2][1]
               xLimLoLo, xLimLoUp = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               xLimUpLo, xLimUpUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1]), 1.
               local evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoLo,xLimLoUp}, idxS[2]==numCells[2]) end
               local evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimUpLo,xLimUpUp}, idxS[2]==numCells[2]) end
               yLimLo = -1.
               yLimUp = wrapNum(cellLower[2] + yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
               -- Add scenario sii-like contribution.
               yLimLo = yLimUp    -- Keep the order of yLimUp and yLimLo here.
               yLimUp = 1.
               y0Lo = cellBounds[2][2]
               xLimLo = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
               xLimUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
               dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
               projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
               projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
               interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
            end

         end

      elseif nilNum == 2 then
         -- There are 6 possible scenarios:
         --   six)   yLowerDest+yShift crosses lower and upper y-boundaries (increasing yShift).
         --   sx)    yLowerDest+yShift crosses lower and upper y-boundaries (decreasing yShift).
         --   sxi)   yUpperDest+yShift crosses lower and upper y-boundaries (decreasing yShift).
         --   sxii)  yUpperDest+yShift crosses lower and upper y-boundaries (increasing yShift).
         --   sxiii) yLowerDest+yShift crosses lower y-boundary and yUpperDest+yShift crosses upper y-boundary (decreasing yShift).
         --   sxiv)  yLowerDest+yShift crosses lower y-boundary and yUpperDest+yShift crosses upper y-boundary (increasing yShift).

         -- The first four of these we will do with a single integral with fixed y-limits, and
         -- a fixed x-limit and a variable x-limit. The last 2 might be easier to do by integrating
         -- the full cell, and subtracting the subregions that are not meant to contribute.

         
         local nilNumDestLo, nilIdxsDestLo, nonNilIdxsDestLo = nilHowmany2(trapCornerX[1])
         local nilNumDestUp, nilIdxsDestUp, nonNilIdxsDestUp = nilHowmany2(trapCornerX[2])
         if nilNumDestLo==2 or nilNumDestUp==2 then
            -- Scenario six, sx, sxi or sxii.
            local yLimLo, yLimUp = -1., 1.
            local evaluateLo, evaluateUp
            if nilNumDestLo==0 then
               -- Scenario six or sx.
               local y0Lo = cellBounds[2][1]
               if trapCornerX[1][1] < trapCornerX[1][2] then
                  -- Scenario six (similar to scenario si).
                  xLimLo, xLimUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1])
                  evaluateLo = function(tCurr,xn) return -1.0 end
                  evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               else
                  -- Scenario sx.
                  xLimLo, xLimUp = (trapCornerX[1][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
                  evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
                  evaluateUp = function(tCurr,xn) return 1.0 end
               end
            else
               local y0Up = cellBounds[2][2]
               if trapCornerX[2][2] < trapCornerX[2][1] then
                  -- Scenario sxi.
                  xLimLo, xLimUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1])
                  evaluateLo = function(tCurr,xn) return -1.0 end
                  evaluateUp = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               else
                  -- Scenario sxii.
                  xLimLo, xLimUp = (trapCornerX[2][1] - xcS[1])/(0.5*dxS[1]), (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
                  evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
                  evaluateUp = function(tCurr,xn) return 1.0 end
               end
            end
            dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
            projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
            projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
            interpGenSub(qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
         else
            -- Scenario sxiii or sxiv.
            local evaluateLo, evaluateUp, evaluateLoP, evaluateUpP
            local yLimUp, yLimLo, yLimUpP, yLimLoP
--            if trapCornerX[2][1]==nil then
--            if yShiftFunc(srcCellBounds[1][1]) > yShiftFunc(trapCornerX[nonNilIdxs[1][1]][nonNilIdxs[1][2]]) then
            if yShiftFunc(srcCellBounds[1][1]) > yShiftFunc(srcCellBounds[1][2]) then
               -- Scenario sxiii (kernel takes information of scenario sii and siii-like integrals).
               -- Scenario sii-like part.
               local y0Lo = cellBounds[2][2]
               yLimLo = wrapNum(cellUpper[2]+yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLo = (yLimLo - xcS[2])/(0.5*dxS[2])
               yLimUp = 1.
               local xLimLo, xLimUp = (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1]), 1.
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
               -- Scenario siii-like part.
               local y0Up = cellBounds[2][1]
               yLimLoP = -1.
               yLimUpP = wrapNum(cellLower[2]+yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUpP = (yLimUpP - xcS[2])/(0.5*dxS[2])
               local xLimLoP, xLimUpP = -1., (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1])
               evaluateLoP = function(tCurr,xn) return -1.0 end
               evaluateUpP = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLoP,xLimUpP}, idxS[2]==numCells[2]) end

            else
               -- Scenario sxiv (kernel takes information of scenario si and siv-like integrals).
               -- Scenario si-like part.
               local y0Lo = cellBounds[2][2]
               yLimLoP = wrapNum(cellUpper[2]+yShiftFunc(srcCellBounds[1][1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimLoP = (yLimLoP - xcS[2])/(0.5*dxS[2])
               yLimUpP = 1.
               local xLimLoP, xLimUpP = -1., (trapCornerX[2][2] - xcS[1])/(0.5*dxS[1])
               evaluateLoP = function(tCurr,xn) return -1.0 end
               evaluateUpP = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Lo, xcS, xc, dx, {xLimLoP,xLimUpP}, idxS[2]==numCells[2]) end
               -- Scenario siv-like part.
               local y0Up = cellBounds[2][1]
               yLimLo = -1.
               yLimUp = wrapNum(cellLower[2]+yShiftFunc(srcCellBounds[1][2]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
               yLimUp = (yLimUp - xcS[2])/(0.5*dxS[2])
               xLimLo, xLimUp = (trapCornerX[1][1] - xcS[1])/(0.5*dxS[1]), 1.
               evaluateLo = function(tCurr,xn) return yShiftNormInvFuncPartialy(xn[1], y0Up, xcS, xc, dx, {xLimLo,xLimUp}, idxS[2]==numCells[2]) end
               evaluateUp = function(tCurr,xn) return 1.0 end
            end
            dyPartial[1], ycPartial[1] = yLimUp-yLimLo, 0.5*(yLimUp+yLimLo)
            projXlimGen(evaluateLo, dyPartial, ycPartial, qLo_y)
            projXlimGen(evaluateUp, dyPartial, ycPartial, qUp_y)
            dyPartialP[1], ycPartialP[1] = yLimUpP-yLimLoP, 0.5*(yLimUpP+yLimLoP)
            projXlimGen(evaluateLoP, dyPartialP, ycPartialP, qLoP_y)
            projXlimGen(evaluateUpP, dyPartialP, ycPartialP, qUpP_y)
            interpGenSubMtwoCorners(qLoP_y:data(), qUpP_y:data(), yLimLoP, yLimUpP, dyPartialP[1], ycPartialP[1], qLo_y:data(), qUp_y:data(), yLimLo, yLimUp, dyPartial[1], ycPartial[1], srcPtr:data(), destPtr:data())
         end

      elseif nilNum == 4 then
         -- There are 2 possible scenarios:
         --   sxv)  yLowerDest+yShift crosses lower and upper x-boundaries.
         --   sxvi) yUpperDest+yShift crosses lower and upper x-boundaries.

         local xLimLo, xLimUp = -1., 1.
         local evaluateLo, evaluateUp
         -- These integrals will use fixed x-limits and variable y limits.
         local yShiftLoXc = wrapNum(cellLower[2]+yShiftFunc(xcS[1]),grid:lower(2),grid:upper(2),idxS[2]==numCells[2])
         if (srcCellBounds[2][1] <= yShiftLoXc) and (yShiftLoXc <= srcCellBounds[2][2]) then
            -- Scenario sxv.
            local y0Lo = cellBounds[2][1]
            evaluateLo = function(tCurr,xn) return yShiftNormFunc(xn[1],y0Lo,xcS,xc,dx,idxS[2]==numCells[2]) end
            evaluateUp = function(tCurr,xn) return 1.0 end
         else
            -- Scenario sxvi.
            local y0Up = cellBounds[2][2]
            evaluateLo = function(tCurr,xn) return -1.0 end
            evaluateUp = function(tCurr,xn) return yShiftNormFunc(xn[1],y0Up,xcS,xc,dx,idxS[2]==numCells[2]) end
         end
         dxPartial[1], xcPartial[1] = xLimUp-xLimLo, 0.5*(xLimUp+xLimLo)
         projXlimGen(evaluateLo, dxPartial, xcPartial, qLo_x)
         projXlimGen(evaluateUp, dxPartial, xcPartial, qUp_x)
         interpGenSubY(xLimLo, xLimUp, qLo_x:data(), qUp_x:data(), dxPartial[1], xcPartial[1], srcPtr:data(), destPtr:data())


      end
--      end

   end
--   end
--   end

end
end

fldDest:write("fldDest.bp")

local intFldDest = DataStruct.DynVector {
   numComponents = 1,
}
intQuant:advance(0., {fldDest}, {intFldDest})
intFldDest:write("intFldDest.bp", 0., 0)
