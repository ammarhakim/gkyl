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
local SerendipityNodes     = require "Lib.SerendipityNodes"
local ffi                  = require "ffi"
local ffiC                 = ffi.C
local xsys                 = require "xsys"
local TwistShiftInterpDecl = require "Updater.interpolateCalcData.TwistShiftInterpModDecl"
GKYL_EMBED_INP             = false
local math                 = require "sci.math"  -- For sign function.

-- Create two fields on a 2D grid. Then interpolate one field onto the other
-- but with a shift in y that is a function of x (assuming periodicity in y).

local polyOrder    = 1
local lower        = {-2.0, -1.50}
local upper        = { 2.0,  1.50}
local numCells     = {20, 20}
local periodicDirs = {2}

local grid = Grid.RectCart {
   lower        = lower,
   upper        = upper,
   cells        = numCells,
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
local fldSrc = DataStruct.Field {
   onGrid        = grid,
   numComponents = basis:numBasis(),
   ghost         = {1, 1},
   metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
}
local fldSrcShifted = DataStruct.Field {
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
   lower        = {lower[1]},
   upper        = {upper[1]},
   cells        = {numCells[1]},
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
                      return -0.3*x+1.
--                      return 0.1*(x+2.)^2+0.2*x+0.6
                   end
project1D:setFunc(function(t, xn) return yShiftFunc(xn[1]) end)
project1D:advance(0., {}, {yShift})
yShift:write("yShift.bp")

local wrapNum = function (val, lower, upper, pickUpper)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   --local pickUpper = pickUpper==nil and false or pickUpper
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


local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn) return 1. end,
}
-- Donor field.
local fldDoFunc = function(t, xn)
   local x, y       = xn[1], xn[2]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))
--   return math.sin((2.*math.pi/3.)*y)
--   return math.exp(-((y-muY)^2)/(2.*(sigY^2)))
--   return 1.
--   if y < 0. then
--      return y+2
--   else
--      return -y+2
--   end
end
-- Shifted donor field.
local fldDoShiftedFunc = function(t, xn)
   local x, y       = xn[1], xn[2]
   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((wrapNum(y+yShiftFunc(x),grid:lower(2),grid:upper(2),true)-muY)^2)/(2.*(sigY^2)))
--   return math.sin((2.*math.pi/3.)*((wrapNum(y+yShiftFunc(x),grid:lower(2),grid:upper(2),true))))
--   return math.exp(-((wrapNum(y+yShiftFunc(x),grid:lower(2),grid:upper(2),true)-muY)^2)/(2.*(sigY^2)))
--   return 1.
--   if y < -1. or y > 0.5 then
--      return y+2
--   else
--      return -y+2
--   end
end

-- Project donor field function onto basis.
project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
project:advance(0., {}, {fldSrc})
fldSrc:write("fldSrc.bp")
-- Project shifted donor field function onto basis.
project:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
project:advance(0., {}, {fldSrcShifted})
fldSrcShifted:write("fldSrcShifted.bp")

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
   assert( xRoot==nil, string.format("yShiftFunc is zero at %f. Please select a function that is not zero in the domain.",xRoot) )
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
-- If y+yShift-ySrc=0 has no roots, it is possible that y+yShift intersects a periodic
-- copy of this domain. Check for such cases by looking for the roots of
-- y+yShift-(ySrc+N*Ly)=0 where Ly is the length of the domain along y and N is an integer.
local yShiftedyIntersect = function(xBounds, yDest, ySrc, jDest, jSrc)
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
   yS = wrapNum(yS,grid:lower(2),grid:upper(2),pickUpper)
   local ySN   = (yS-xcSrc[2])/(0.5*dxDest[2])
   return ySN
end

-- Given a normalized yShift, the source and destination cell centers and the
-- destination cell length, compute the normalized x (i.e. in [-1,1]).
local yShiftNormInvFunc = function(ySNIn, y0In, xcSrcIn, xcDestIn, dxDestIn, xLim, pickUpper)
   -- y+yShift intersects the x-boundary instead of the y-boundary. This means that the y+yShift
   -- line that the defines the x-intergral limit is contained within a subset of the y-range
   -- of the source cell.
   -- This function inverts that y+yShifted(x) line. This function gets passed to a projection
   -- operation, so that the resulting DG expansion is defined in a subcell of the source cell.
   local function lossF(xNIn)
      return ySNIn - yShiftNormFunc(xNIn, y0In, xcSrcIn, xcDestIn, dxDestIn, pickUpper)
   end
   local tol = 1.e-11
   -- Have to consider extended logical space even though the part of the y+yShift curve
   -- that is not in this cell does not get used in the integral because the fixed
   -- y limit cuts if off. But if we didn't consider extended logical space the shape of the
   -- x-limit function would not come out right.
   local eps = 0.  -- Need wiggle room because it might fail if the root is at the boundary.
   local xOutN = root.ridders(lossF, xLim[1]-eps, xLim[2]+eps, rootStop(tol))
   return xOutN
end

local doTarOff = function(dxSrc, xcSrc, xcDest)
   -- y-offset between the donor and the target cell, in the direction of the shift.
   local shift     = yShiftFunc(xcSrc[1])
   local shiftSign = sign(shift)
   local Ly        = grid:upper(2)-grid:lower(2)
   local yShLyFac  = math.floor(math.abs(shift)/Ly)

   if yShLyFac < 1. then yShLyFac = 1. end

   local ps = 0.
   if shiftSign==1 and xcSrc[2] <= xcDest[2] then
      return xcSrc[2] + yShLyFac*Ly - xcDest[2] + ps
   elseif shiftSign==-1 and xcSrc[2] >= xcDest[2] then
      return xcSrc[2] - yShLyFac*Ly - xcDest[2] + ps
   else
      return xcSrc[2] - xcDest[2] + ps
   end
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

local p2l = function(valIn, xcIn, dxIn)
   -- Transform a physical coordinate (valIn) to the [-1,1] logical
   -- space in a cell with cell center xcIn and length dxIn.
   return 2.*(valIn-xcIn)/dxIn
end
local dxAxc = function(aIn,bIn)
   -- Return the length and cell center of the interval [a,b].
   return bIn-aIn, 0.5*(aIn+bIn)
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
local etaLo_xi, etaUp_xi   = Lin.Vec(numNodes), Lin.Vec(numNodes)
local xiLo_eta, xiUp_eta   = Lin.Vec(numNodes), Lin.Vec(numNodes)
local xiLoL_eta, xiUpL_eta = Lin.Vec(numNodes), Lin.Vec(numNodes)
local xiLoR_eta, xiUpR_eta = Lin.Vec(numNodes), Lin.Vec(numNodes)
local compToPhys = loadstring(compToPhysTempl {NDIM = 1} )()
local xi = Lin.Vec(1)                  -- Coordinates at node.
local dxSub, xcSub   = Lin.Vec(1), Lin.Vec(1)
local dySub, ycSub   = Lin.Vec(1), Lin.Vec(1)
local dySubL, ycSubL = Lin.Vec(1), Lin.Vec(1)
local dySubR, ycSubR = Lin.Vec(1), Lin.Vec(1)
local projLim = function(funcIn, dxIn, xcIn, q_yOut)
   for i = 1, numNodes do
      compToPhys(nodes[i], dxIn, xcIn, xi)      -- Compute physical coordinate xi.
      evalFunc(0.0, xi, funcIn, fv[i])
   end
   ffiC.nodToMod(fv:data(), numNodes, numVal, 1, polyOrder1D, q_yOut:data())
end

local intSubX           = TwistShiftInterpDecl.selectTwistShiftInterpSubX(basis:id(), 2, polyOrder1D)
local intSubY           = TwistShiftInterpDecl.selectTwistShiftInterpSubY(basis:id(), 2, polyOrder1D)
local intSubMtwoCorners = TwistShiftInterpDecl.selectTwistShiftInterpSubMinusTwoCorners(basis:id(), 2, polyOrder1D)

local subCellInt_xLimDG = function(xLimFuncs,yBounds,dyDo,offDoTar)
   -- Perform a sub-cell integral with x-limits (possibly functions of y) given by a
   -- DG polynomial representation. Here x-y mean the logical x-y of the donor cell.
   -- The functions xLoLim/xUpLim define these x-limits in the y-logical space of the
   -- donor cell. Then we project these functions onto a 1D DG basis. Since the integral
   -- may not span the whole y-logical space of the donor cell (dySub<2), the 1D
   -- projection is done treating the partial y-logical space of the donor cell as a
   -- 'physical' space on which the x-limit polynomial is defined, which a different
   -- relationship to its own logical coordinate than that between the logical and physical
   -- coordinate of the donor cell. This means that the difference between the logical
   -- coordinates of the x-limits and the donor field needs to be accounted for. We do that
   -- using dySub and ycSub in the kernel.
   --   xLimFuncs: functions defining lower/upper x-limits of integral.
   --   yBounds:   lower/upper logical y-limits of the integral.
   --   dy:        physical y-cell length of the donor cell.
   --   offDoTar:  physical y-offset between donor and target cells.
   dySub[1], ycSub[1] = dxAxc(yBounds[1],yBounds[2])
   if dySub[1] > 0 then 
      projLim(xLimFuncs[1], dySub, ycSub, xiLo_eta)
      projLim(xLimFuncs[2], dySub, ycSub, xiUp_eta)
      intSubX(xiLo_eta:data(), xiUp_eta:data(), yBounds[1], yBounds[2], dySub[1], ycSub[1],
              dyDo, offDoTar, yShiftPtr:data(), srcPtr:data(), destPtr:data())
   end
end

local subcellTrapezoidInt = function(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, pickUpper)
   -- Perform integral over a trapezoidal subcell region.
   --   xLogBounds: 2x2 table with the bounds of the logical x lower and upper limits.
   --   yLogBounds: 2 element table with the logical y-bounds of the logical x lower and upper limits.
   --   yTars: y_{j-/+1/2} in target cell which defines the lower/upper x limits.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.

   local ycOff = doTarOff(dxOr, xcOr, xcTar)
   -- Functions describing the lower and upper limits in logical [-1,1] x-space.
   local xLogLims = {
      function(t,xn)
         return yShiftNormInvFunc(xn[1], yTars[1], xcOr, xcTar, dxTar, xLogBounds[1], pickUpper)
      end,
      function(t,xn)
         return yShiftNormInvFunc(xn[1], yTars[2], xcOr, xcTar, dxTar, xLogBounds[2], pickUpper)
      end
   }
   subCellInt_xLimDG(xLogLims,yLogBounds,dxOr[2],ycOff)
end

local subCellInt_sN = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sN subcell integral.
   -- The subcell region is a 4-sided trapezoid. Both y_{j -/+ 1/2}+yShift
   -- intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   local yLogBounds = {-1., 1.}
   -- Establish if the y_{j-/+1/2}+yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=yLowerSrc crosses y_{j-/+1/2}+yShift.
   local yTars, xLogBounds = {nil,nil}, {{nil,nil},{nil,nil}}
   if x_pq[1][1] < x_pq[2][1] then   -- Monotonically decreasing yShift.
      yTars[1], yTars[2] = boundsTar[2][1], boundsTar[2][2]
      xLogBounds[1][1], xLogBounds[1][2] = p2l(x_pq[1][2],xcOr[1],dxOr[1]), p2l(x_pq[1][1],xcOr[1],dxOr[1])
      xLogBounds[2][1], xLogBounds[2][2] = p2l(x_pq[2][2],xcOr[1],dxOr[1]), p2l(x_pq[2][1],xcOr[1],dxOr[1])
   else                              -- Monotonically increasing yShift.
      yTars[1], yTars[2] = boundsTar[2][2], boundsTar[2][1]
      xLogBounds[1][1], xLogBounds[1][2] = p2l(x_pq[2][1],xcOr[1],dxOr[1]), p2l(x_pq[2][2],xcOr[1],dxOr[1])
      xLogBounds[2][1], xLogBounds[2][2] = p2l(x_pq[1][1],xcOr[1],dxOr[1]), p2l(x_pq[1][2],xcOr[1],dxOr[1])
   end
   subcellTrapezoidInt(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
end

local subCellInt_xUpLimDG = function(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, pickUpper)
   -- Perform subcell integral with variable upper x-limit.
   --   yTar: yTar+yShift gives the curve defining the upper x-limit.
   --   xLogBounds: logical x range in which to invert yTar+yShift.
   --   yLogBounds: logical y limits of the integral.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   local ycOff    = doTarOff(dxOr, xcOr, xcTar)
   local xLogLims = {function(t,xn) return -1.0 end,
                     function(t,xn)
                        return yShiftNormInvFunc(xn[1], yTar, xcOr, xcTar, dxTar, xLogBounds, pickUpper)
                     end}
   subCellInt_xLimDG(xLogLims,yLogBounds,dxOr[2],ycOff)
end

local subCellInt_xLoLimDG = function(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, pickUpper)
   -- Perform subcell integral with variable lower x-limit.
   --   yTar: yTar+yShift gives the curve defining the upper x-limit.
   --   xLogBounds: logical x range in which to invert yTar+yShift.
   --   yLogBounds: logical y limits of the integral.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   local ycOff    = doTarOff(dxOr, xcOr, xcTar)
   local xLogLims = {function(t,xn)
                        return yShiftNormInvFunc(xn[1], yTar, xcOr, xcTar, dxTar, xLogBounds, pickUpper)
                     end,
                     function(t,xn) return 1.0 end}
   subCellInt_xLimDG(xLogLims,yLogBounds,dxOr[2],ycOff)
end

local subCellInt_siORsii = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario si or sii subcell integral.
   -- The subcell region is 3-sided, abutting one of the upper corners of the source cell.
   -- Only y_{j - 1/2}+yShift intersects y_{j+m + 1/2} (the upper y-boundary of srcCell).
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each dierection
   local yShiftUp = yShiftFunc(boundsTar[1][2])
   local yTar     = boundsTar[2][1]
   if yShiftUp >= yShiftFunc(x_pq[1][2]) then
      -- Scenario si.
      local yLogBounds = {p2l(wrapNum(boundsTar[2][1]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2]), 1.}
      local xLogBounds = {-1, p2l(x_pq[1][2],xcOr[1],dxOr[1])}
      subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
   else
      -- Scenario sii.
      local yLogBounds = {p2l(wrapNum(boundsTar[2][1]+yShiftUp,grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2]), 1.}
      local xLogBounds = {p2l(x_pq[1][2],xcOr[1],dxOr[1]), 1}
      subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
   end
end

local subCellInt_siiiORsiv = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario siii or siv subcell integral.
   -- The subcell region is 3-sided, abutting one of the lower corners of the source cell.
   -- Only y_{j + 1/2}+yShift intersects y_{j+m - 1/2} (the upper y-boundary of srcCell).
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each dierection
   local yShiftUp = yShiftFunc(boundsTar[1][2])
   local yTar     = boundsTar[2][2]
   if yShiftUp <= yShiftFunc(x_pq[2][1]) then
      -- Scenario siii.
      local yLogBounds = {-1., p2l(wrapNum(boundsTar[2][2]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])}
      local xLogBounds = {-1, p2l(x_pq[2][1],xcOr[1],dxOr[1])}
      subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
   else
      -- Scenario siv.
      local yLogBounds = {-1., p2l(wrapNum(boundsTar[2][2]+yShiftUp,grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])}
      local xLogBounds = {p2l(x_pq[2][1],xcOr[1],dxOr[1]), 1}
      subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
   end
end

local subCellInt_six = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario six subcell integral.
   --   x_pq: points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr: cell center coordinates of origin cell.
   --   dxOr: cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells: number of cells in each dierection
   local xLogBounds = {p2l(x_pq[1][1],xcOr[1],dxOr[1]), p2l(x_pq[1][2],xcOr[1],dxOr[1])}
   local yLogBounds = {-1., 1.}
   local yTar       = boundsTar[2][1]    -- y_{j -/+ 1/2} to shift in y.
   subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
end

local subCellInt_sx = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sx subcell integral.
   --   x_pq: points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr: cell center coordinates of origin cell.
   --   dxOr: cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells: number of cells in each dierection
   local xLogBounds = {p2l(x_pq[1][2],xcOr[1],dxOr[1]), p2l(x_pq[1][1],xcOr[1],dxOr[1])}
   local yLogBounds = {-1., 1.}
   local yTar       = boundsTar[2][1]    -- y_{j -/+ 1/2} to shift in y.
   subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
end

local subCellInt_sxi = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sxi subcell integral.
   --   x_pq: points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr: cell center coordinates of origin cell.
   --   dxOr: cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells: number of cells in each dierection
   local xLogBounds = {p2l(x_pq[2][2],xcOr[1],dxOr[1]), p2l(x_pq[2][1],xcOr[1],dxOr[1])}
   local yLogBounds = {-1., 1.}
   local yTar       = boundsTar[2][2]    -- y_{j -/+ 1/2} to shift in y.
   subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
end

local subCellInt_sxii = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sxii subcell integral.
   --   x_pq: points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr: cell center coordinates of origin cell.
   --   dxOr: cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells: number of cells in each dierection
   local xLogBounds = {p2l(x_pq[2][1],xcOr[1],dxOr[1]), p2l(x_pq[2][2],xcOr[1],dxOr[1])}
   local yLogBounds = {-1., 1.}
   local yTar       = boundsTar[2][2]
   subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
end

local subCellInt_svORsvi = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sv or svi subcell integral.
   -- The subcell region is 5-sided, abutting one of the lower corners of the source cell.
   -- Only y_{j + 1/2}+yShift does not intersects y_{j+m - 1/2} (the lower y-boundary of srcCell).
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each dierection
   -- Establish if the y_{j-/+1/2}+yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=yUpperSrc crosses y_{j+m-/+1/2}+yShift.
   if x_pq[1][2] > x_pq[2][2] then   -- Scenario sv.
      local xLogLoLo, xLogLoUp = -1., p2l(x_pq[2][2],xcOr[1],dxOr[1])
      if (xLogLoUp-xLogLoLo) > 0 then
         -- Add scenario sN-like contribution.
         local yTars      = {boundsTar[2][2], boundsTar[2][1]}
         local xLogBounds = {{xLogLoLo, xLogLoUp},
                             {p2l(x_pq[1][1],xcOr[1],dxOr[1]), p2l(x_pq[1][2],xcOr[1],dxOr[1])}}
         local yLogBounds = {p2l(wrapNum(boundsTar[2][2]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2]), 1.}
         subcellTrapezoidInt(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
         -- Add scenario siii-like contribution.
         local yTar    = boundsTar[2][1]
         xLogBounds    = {p2l(x_pq[1][1],xcOr[1],dxOr[1]), p2l(x_pq[1][2],xcOr[1],dxOr[1])}
         yLogBounds[2] = yLogBounds[1]    -- Keep the order of yLogBounds here.
         yLogBounds[1] = -1.
         subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
      else   -- This is better handled as a scenario six.
         subCellInt_six(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
      end
   else   -- Scenario svi.
      local xLogUpLo, xLogUpUp = p2l(x_pq[2][2],xcOr[1],dxOr[1]), 1.
      if (xLogUpUp-xLogUpLo) > 0 then
         -- Add scenario sN-like contribution.
         local yTars      = {boundsTar[2][1], boundsTar[2][2]} 
         local xLogBounds = {{p2l(x_pq[1][2],xcOr[1],dxOr[1]), p2l(x_pq[1][1],xcOr[1],dxOr[1])},
                             {xLogUpLo, xLogUpUp}}
         local yLogBounds = {p2l(wrapNum(boundsTar[2][2]+yShiftFunc(boundsTar[1][2]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2]), 1.} 
         subcellTrapezoidInt(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
         -- Add scenario siv-like contribution.
         local yTar    = boundsTar[2][1]
         xLogBounds    = {p2l(x_pq[1][2],xcOr[1],dxOr[1]), p2l(x_pq[1][1],xcOr[1],dxOr[1])}
         yLogBounds[2] = yLogBounds[1]    -- Keep the order of yLogBounds here.
         yLogBounds[1] = -1. 
         subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
      else    -- This is better handled as a scenario sx.
         subCellInt_sx(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
      end
   end
end

local subCellInt_sviiORsviii = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario svii or sviii subcell integral.
   -- The subcell region is 5-sided, abutting one of the upper corners of the source cell.
   -- Only y_{j - 1/2}+yShift does not intersects y_{j+m + 1/2} (the upper y-boundary of srcCell).
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each direction
   -- Establish if the y_{j-/+1/2}+yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=yLowerSrc crosses y_{j+m-/+1/2}+yShift.
   if x_pq[1][1] < x_pq[2][1] then   -- Scenario svii.
      local xLimLoLo, xLimLoUp = -1., p2l(x_pq[1][1],xcOr[1],dxOr[1])
      if (xLimLoUp-xLimLoLo) > 0 then
         -- Add scenario sN-like contribution.
         local yTars      = {boundsTar[2][1], boundsTar[2][2]} 
         local xLogBounds = {{xLimLoLo, xLimLoUp},
                             {p2l(x_pq[2][2],xcOr[1],dxOr[1]), p2l(x_pq[2][1],xcOr[1],dxOr[1])}}
         local yLogBounds = {-1., p2l(wrapNum(boundsTar[2][1]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])} 
         subcellTrapezoidInt(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
         -- Add scenario si-like contribution.
         local yTar    = boundsTar[2][2]
         xLogBounds    = {p2l(x_pq[2][2],xcOr[1],dxOr[1]), p2l(x_pq[2][1],xcOr[1],dxOr[1])}
         yLogBounds[1] = yLogBounds[2]    -- Keep the order of yLogBounds here.
         yLogBounds[2] = 1.
         subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
      else   -- This is better handled as a scenario sxi.
         subCellInt_sxi(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
      end
   else   -- Scenario sviii.
      local xLimUpLo, xLimUpUp = p2l(x_pq[1][1],xcOr[1],dxOr[1]), 1.
      if (xLimUpUp-xLimUpLo) > 0 then
         -- Add scenario sN-like contribution.
         local yTars      = {boundsTar[2][2], boundsTar[2][1]}
         local xLogBounds = {{p2l(x_pq[2][1],xcOr[1],dxOr[1]), p2l(x_pq[2][2],xcOr[1],dxOr[1])},
                             {xLimUpLo, xLimUpUp}}
         local yLogBounds = {-1., p2l(wrapNum(boundsTar[2][1]+yShiftFunc(boundsTar[1][2]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])}
         subcellTrapezoidInt(xLogBounds, yLogBounds, yTars, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
         -- Add scenario sii-like contribution.
         local yTar       = boundsTar[2][2]
         local xLogBounds = {p2l(x_pq[2][1],xcOr[1],dxOr[1]), p2l(x_pq[2][2],xcOr[1],dxOr[1])} 
         yLogBounds[1]    = yLogBounds[2]    -- Keep the order of yLogBounds here.
         yLogBounds[2]    = 1.
         subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xcOr, dxOr, xcTar, dxTar, idxOr[2]==numCells[2])
      else   -- This is better handled as a scenario sxii.
         subCellInt_sxii(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
      end
   end
end

local subCellInt_sxiiiORsxiv = function(x_pq, idxOr, xcOr, dxOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario svxiii or sxiv subcell integral.
   -- The subcell region is 5-sided, abutting one of the upper corners of the source cell.
   -- Only y_{j - 1/2}+yShift does not intersects y_{j+m + 1/2} (the upper y-boundary of srcCell).
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each direction
   local ycOff = doTarOff(dxOr, xcOr, xcTar)
   local yTarL, yTarR
   local xLogLimsL, xLogLimsR     = {}, {}
   local xLogBoundsL, xLogBoundsR = {}, {}
   local yLogBoundsL, yLogBoundsR = {}, {}
   if yShiftFunc(boundsTar[1][1]) < yShiftFunc(boundsTar[1][2]) then   -- Scenario sxiii.
      -- Scenario si-like part.
      yTarL       = boundsTar[2][2]
      xLogBoundsL = {-1., p2l(x_pq[2][2],xcOr[1],dxOr[1])}
      yLogBoundsL = {p2l(wrapNum(boundsTar[2][2]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2]), 1.}
      -- Scenario siv-like part.
      yTarR       = boundsTar[2][1]
      xLogBoundsR = {p2l(x_pq[1][1],xcOr[1],dxOr[1]), 1.}
      yLogBoundsR = {-1., p2l(wrapNum(boundsTar[2][1]+yShiftFunc(boundsTar[1][2]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])}
   else   -- Scenario sxiv.
      -- Scenario sii-like part.
      yTarR       = boundsTar[2][2]
      xLogBoundsR = {p2l(x_pq[2][2],xcOr[1],dxOr[1]), 1.}
      yLogBoundsR = {p2l(wrapNum(boundsTar[2][2]+yShiftFunc(boundsTar[1][2]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]), xcOr[2], dxOr[2]), 1.}
      -- Scenario siii-like part.
      yTarL       = boundsTar[2][1]
      xLogBoundsL = {-1., p2l(x_pq[1][1],xcOr[1],dxOr[1])}
      yLogBoundsL = {-1., p2l(wrapNum(boundsTar[2][1]+yShiftFunc(boundsTar[1][1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2]),xcOr[2],dxOr[2])}
   end
   -- Project x-limits of left subcell region to subtract.
   xLogLimsL = {function(t,xn) return -1.0 end,
                function(t,xn)
                   return yShiftNormInvFunc(xn[1], yTarL, xcOr, xcTar, dxTar, xLogBoundsL, idxOr[2]==numCells[2])
                end}
   dySubL[1], ycSubL[1] = dxAxc(yLogBoundsL[1],yLogBoundsL[2])
   projLim(xLogLimsL[1], dySubL, ycSubL, xiLoL_eta)
   if (xLogBoundsL[2]-xLogBoundsL[1]) > 1.e-14 and dySubL[1] > 1.e-14 then
      projLim(xLogLimsL[2], dySubL, ycSubL, xiUpL_eta)
   else
      projLim(xLogLimsL[1], dySubL, ycSubL, xiUpL_eta)
   end
   -- Project x-limits of right subcell region to subtract.
   xLogLimsR = {function(t,xn)
                   return yShiftNormInvFunc(xn[1], yTarR, xcOr, xcTar, dxTar, xLogBoundsR, idxOr[2]==numCells[2])
                end,
                function(t,xn) return 1.0 end}
   dySubR[1], ycSubR[1] = dxAxc(yLogBoundsR[1],yLogBoundsR[2])
   if (xLogBoundsR[2]-xLogBoundsR[1]) > 0 and dySubR[1] > 1.e-14  then
      projLim(xLogLimsR[1], dySubR, ycSubR, xiLoR_eta)
   else
      projLim(xLogLimsR[2], dySubR, ycSubR, xiLoR_eta)
   end
   projLim(xLogLimsR[2], dySubR, ycSubR, xiUpR_eta)
   -- The MtwoCorners kernel divides by dy, so we need to set them to something other than 0
   -- to avoid NaNs. The corresponding contribution is zero anyway because the integral limits are equal.
   if dySubL[1] < 1.e-14 then dySubL[1] = 1. end
   if dySubR[1] < 1.e-14 then dySubR[1] = 1. end
   intSubMtwoCorners(xiLoL_eta:data(), xiUpL_eta:data(), yLogBoundsL[1], yLogBoundsL[2], dySubL[1], ycSubL[1],
                     xiLoR_eta:data(), xiUpR_eta:data(), yLogBoundsR[1], yLogBoundsR[2], dySubR[1], ycSubR[1],
                     dxOr[2], ycOff, yShiftPtr:data(), srcPtr:data(), destPtr:data())
end

local subCellInt_sxvORsxvi = function(x_pq, idxOr, xcOr, dxOr, boundsOr, xcTar, dxTar, boundsTar, numCells)
   -- Perform scenario sxv or sxvi subcell integral.
   -- These integrals will use fixed x-limits and variable y limits.
   --   x_pq:  points where y_{j -/+ 1/2}+yShift intersect y_{j+m -/+ 1/2} (the upper/lower y-boundaries of srcCell).
   --   idxOr: multi-dimensional index of origin cell.
   --   xcOr:  cell center coordinates of origin cell.
   --   dxOr:  cell lengths of origin cell.
   --   xcTar: cell center coordinates of target cell.
   --   dxTar: cell lengths of target cell.
   --   boundsTar: cell boundaries along each direction of target cell.
   --   numCells:  number of cells in each direction
   local ycOff      = doTarOff(dxOr, xcOr, xcTar)
   local xLogBounds = {-1., 1.}
   local xLogLims   = {}
   local yShiftLoXc = wrapNum(boundsTar[2][1]+yShiftFunc(xcOr[1]),grid:lower(2),grid:upper(2),idxOr[2]==numCells[2])
   if (boundsOr[2][1] <= yShiftLoXc) and (yShiftLoXc <= boundsOr[2][2]) then   -- Scenario sxv.
      local yTar = boundsTar[2][1]
      xLogLims = {function(t,xn)
                     return yShiftNormFunc(xn[1],yTar,xcOr,xcTar,dxTar,idxOr[2]==numCells[2])
                  end,
                  function(t,xn) return 1.0 end}
   else   -- Scenario sxvi.
      local yTar = boundsTar[2][2]
      xLogLims = {function(t,xn) return -1.0 end,
                  function(t,xn)
                     return yShiftNormFunc(xn[1],yTar,xcOr,xcTar,dxTar,idxOr[2]==numCells[2])
                  end}
   end
   dxSub[1], xcSub[1] = dxAxc(xLogBounds[1],xLogBounds[2])
   projLim(xLogLims[1], dxSub, xcSub, etaLo_xi)
   projLim(xLogLims[2], dxSub, xcSub, etaUp_xi)
   intSubY(xLogBounds[1], xLogBounds[2], etaLo_xi:data(), etaUp_xi:data(), dxSub[1], xcSub[1],
           dxOr[2], ycOff, yShiftPtr:data(), srcPtr:data(), destPtr:data())
end


local t1 = os.clock()
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
   local delta    = {1.e-9*grid:dx(1),1.e-9*grid:dx(2)}
   for dC = 1,2 do            -- Loop over x=const and y=const boundaries.
      for xS = 1,2 do         -- Loop over lower/upper boundary.

         -- Do the first point separately (it requires pickLower=false in findCell).
         local evPoint = {cellLower[1]+delta[1], cellLower[2]+delta[2]}
         evPoint[dC]   = evPoint[dC]+(xS-1)*(grid:dx(dC)-2.*delta[dC])

         basis1D:evalBasis({(evPoint[1]-grid:cellCenterInDir(1))/(0.5*grid:dx(1))}, basis1Dev)
         local yShiftB = 0.   -- y-shift at boundary.
         for k = 1, numBasis1D do
            yShiftB = yShiftB + yShiftPtr[k]*basis1Dev[k]
         end

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
            local searchPoint = {newP[1],wrapNum(newP[2]+yShiftB,grid:lower(2),grid:upper(2),true)}
            grid:findCell(searchPoint,idxShifted,true,{idx[1],nil})
            if lume.findTable(srcCells,idxShifted)==nil then table.insert(srcCells,idxShifted) end
         end
      end
   end

   if idx[2]==1 then print(string.format("idx = (%d,%d)",idx[1],idx[2])) end

   for iC = 1, #srcCells do
      -- In each contributing cell, approximate the functions qInvL(y) and qInvU(y)
      -- for the min/max of the x-integral with a polynomial. Compute the coefficients
      -- of that polynomial with a projection of yShiftNormInvFunc onto the local basis.
      local idxS = srcCells[iC]
      if idx[2]==1 then print(string.format("  from = (%d,%d)",idxS[1],idxS[2])) end

      grid:setIndex(idxS)
      grid:cellCenter(xcS)
      grid:getDx(dxS)
      fldSrc:fill(indexer(idxS), srcPtr)
      local srcCellBounds = { {grid:cellLowerInDir(1),grid:cellUpperInDir(1)},
                              {grid:cellLowerInDir(2),grid:cellUpperInDir(2)} }

      -- Find the points where y_{j-/+1/2}+yShift intersect the y=y_{j+m-/+1/2} lines.
      local interPts = {{nil, nil}, {nil, nil}}
      for i = 1, 2 do
         for j = 1, 2 do 
            local yDest, ySrc = cellBounds[2][i], srcCellBounds[2][j] 
            interPts[i][j] = yShiftedyIntersect(cellBounds[1], yDest, ySrc, idx[2], idxS[2])
--            if interPts[i][j]==nil then
--               print(string.format("interPts[%d][%d] = ",i,j),nil)
--            else
--               print(string.format("interPts[%d][%d] = %1.12e",i,j,interPts[i][j]))
--            end
         end
      end

      if not nilAny2x2(interPts) then   -- Scenario sN.
         subCellInt_sN(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
      else

         local nilNum, nilIdxs, nonNilIdxs = nilHowmany2x2(interPts)
         
         if nilNum == 3 then
            -- si:   y_{j-1/2}+yShift intersects x_{i-1/2}.
            -- sii:  y_{j-1/2}+yShift intersects x_{i+1/2}.
            -- siii: y_{j+1/2}+yShift intersects x_{i-1/2}.
            -- siv:  y_{j+1/2}+yShift intersects x_{i+1/2}.
            if nonNilIdxs[1][1]==1 then    -- Scenario si or sii.
               subCellInt_siORsii(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
            else   -- Scenario siii or siv.
               subCellInt_siiiORsiv(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
            end
         elseif nilNum == 1 then
            -- sv:    y_{j+1/2}+yShift doesn't intersect y_{j+m-1/2} & intersects x_{i-1/2}.
            -- svi:   y_{j+1/2}+yShift doesn't intersect y_{j+m-1/2} & intersects x_{i+1/2}.
            -- svii:  y_{j-1/2}+yShift doesn't intersect y_{j+m+1/2} & intersects x_{i-1/2}.
            -- sviii: y_{j-1/2}+yShift doesn't intersect y_{j+m+1/2} & intersects x_{i+1/2}.
            if (nilIdxs[1][1]==2) and (nilIdxs[1][2]==1) then   -- Scenario sv or svi.
               subCellInt_svORsvi(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
            else   -- Scenario svii or sviii.
               subCellInt_sviiORsviii(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
            end
         elseif nilNum == 2 then
            -- six:   y_{j-1/2}+yShift crosses y_{j+m-/+1/2} (increasing yShift).
            -- sx:    y_{j-1/2}+yShift crosses y_{j+m-/+1/2} (decreasing yShift).
            -- sxi:   y_{j+1/2}+yShift crosses y_{j+m-/+1/2} (decreasing yShift).
            -- sxii:  y_{j+1/2}+yShift crosses y_{j+m-/+1/2} (increasing yShift).
            -- sxiii: y_{j-1/2}+yShift crosses y_{j+m-1/2} & y_{j+1/2}+yShift crosses y_{j+m+1/2} (increasing yShift).
            -- sxiv:  y_{j-1/2}+yShift crosses y_{j+m-1/2} & y_{j+1/2}+yShift crosses y_{j+m+1/2} (decreasing yShift).
            local nilNumDestLo, nilIdxsDestLo, nonNilIdxsDestLo = nilHowmany2(interPts[1])
            local nilNumDestUp, nilIdxsDestUp, nonNilIdxsDestUp = nilHowmany2(interPts[2])
            if nilNumDestLo==2 or nilNumDestUp==2 then   -- Scenario six, sx, sxi or sxii.
               if nilNumDestLo==0 then   -- Scenario six or sx.
                  if interPts[1][1] < interPts[1][2] then    -- Scenario six (similar to scenario si).
                     subCellInt_six(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
                  else   -- Scenario sx.
                     subCellInt_sx(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
                  end
               else
                  if interPts[2][2] < interPts[2][1] then   -- Scenario sxi.
                     subCellInt_sxi(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
                  else   -- Scenario sxii.
                     subCellInt_sxii(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
                  end
               end
            else   -- Scenario sxiii or sxiv.
               subCellInt_sxiiiORsxiv(interPts, idxS, xcS, dxS, xc, dx, cellBounds, numCells)
            end
         elseif nilNum == 4 then
            -- sxv:  y_{j-1/2}+yShift crosses x_{i-/+1/2}.
            -- sxvi: y_{j+1/2}+yShift crosses x_{i-/+1/2}.
            subCellInt_sxvORsxvi(interPts, idxS, xcS, dxS, srcCellBounds, xc, dx, cellBounds, numCells)
         end
      end

   end

end
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")

fldDest:write("fldDest.bp")

local intFldDest = DataStruct.DynVector {
   numComponents = 1,
}
intQuant:advance(0., {fldDest}, {intFldDest})
intFldDest:write("intFldDest.bp", 0., 0)
