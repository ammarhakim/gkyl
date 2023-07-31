-- Gkyl ------------------------------------------------------------------------
--
-- A module of Lua functions used by the TwistShift updater.
-- 
-- Notes:
--  a) Hard-coded parameters:
--       - wrapNum: eps.
--       - getDonors: stepFrac, deltaFrac.
--       - findIntersect: tol, stepFrac.
--       - yShiftedLogInv: tol, eps.
--  b) At various times we refer to xi and eta, these are logical space coordinates
--     related to physical space coordinates via:
--       x = x_i+(dx/2)*xi
--       y = y_j+(dx/2)*eta
--
-- Module functions:
--   - selectTwistShiftKernels
--   - rootFind(funcIn, lims, tol)
--   - set_yShiftF
--   - set_domLim(gridIn)
--   - set_dx(gridIn)
--   - getDonors
--   - matVec_alloc
--   - set_projData(shiftPolyOrder)
--   - preCalcMat(grid, yShift, doCells)
--
-- Local helper functions.
--   - p2l
--   - dxAxc
--   - wrapNum 
--   - riddersStop
--   - findIntersect
--   - doTarOff
--   - yShiftedLog
--   - yShiftedLogInv
--   - nodToModProj1D
--   - subCellInt_xLimDG
--   - subCellInt_xLoLimDG 
--   - subCellInt_xUpLimDG
--   - subCellInt_yLimDG
--   - subCellInt_trapezoid
--   - subCellInt_sN
--   - subCellInt_siORsii
--   - subCellInt_siiiORsiv
--   - subCellInt_svORsvi
--   - subCellInt_sviiORsviii
--   - subCellInt_sviiORsviii
--   - subCellInt_six
--   - subCellInt_sx
--   - subCellInt_sxi
--   - subCellInt_sxii
--   - subCellInt_sxiiiORsxiv
--   - subCellInt_sxvORsxvi
--   - missingInterPt
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin  = require "Lib.Linalg"
local Mpi  = require "Comm.Mpi"
local lume = require "Lib.lume"
local math = require "sci.math"  -- For sign function.
local root = require "sci.root"
-- The following are needed for projecting onto the basis.
local SerendipityNodes = require "Lib.SerendipityNodes"
local ffi              = require "ffi"
local xsys             = require "xsys"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_bc_twistshift gkyl_bc_twistshift;

/**
 * Create a new updater to apply twist-shift BCs.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded.
 * @param grid Grid dynamic field is defined on.
 * @param cdim Number of configuration space dimensions.
 * @param yshift Discrete y-shift defined on a 1D x-grid.
 * @param ndonors Number of donors for each x-cell.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_new(int dir, int do_dir, int shift_dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const struct gkyl_range *local_range_update, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim,
  const struct gkyl_array *yshift, const int *ndonors, const int *cells_do, bool use_gpu);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the sub-cell integral that
 * has variably x limits represented by a DG polynomial, and a y-integral that
 * goes from yLimLo to yLimUp.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param xLimLo DG representation of the lower x-limit.
 * @param xLimUp DG representation of the upper x-limit.
 * @param yLimLo lower y-limit.
 * @param yLimUp upper y-limit.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
  double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, int cellidx, int doidx);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the sub-cell integral that
 * has variable y limits represented by a DG polynomial, and a x-integral that
 * goes from xLimLo to xLimUp.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param xLimLo lower x-limit.
 * @param xLimUp upper x-limit.
 * @param yLimLo DG representation of the lower y-limit.
 * @param yLimUp DG representation of the upper y-limit.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
  double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, int cellidx, int doidx);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the full-cell integral.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
  double dyDo, double yOff, const double *ySh, int cellidx, int doidx);

/**
 * Multiply the donor matrices by the donor cell dg coefficients
 *
 * @param matsdo nmat of donor matrices for all x grid cells and all donor cells
 * @param vecsdo nmat of donor vectors for all x grid cells and all donor cells
 * @param vecstar nmat of target vectors for all x grid cells
 */
void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Free memory associated with bc_twistshift updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up);

]]






local debugPrint = false  -- Set to false to turn off debugging print statements.
local function printDebug(str)
   if debugPrint then print(str) end
end

-- Kernels. Assigned by the selectTwistShiftKernels function.
local intSubXlimDG = nil 
local intSubYlimDG = nil
local intFullCell  = nil

local _M = {}

function _M.init(dir, do_dir, shift_dir, edge, localExtRange, localUpdateRange, num_ghosts, basis, grid, cdim, yshift, ndonors, cells_do, use_gpu)
   _M._zero = ffi.gc(ffiC.gkyl_bc_twistshift_new(dir, do_dir, shift_dir, edge, localExtRange, localUpdateRange, num_ghosts:data(), basis._zero, grid._zero, cdim, yshift, ndonors:data(), cells_do:data(), use_gpu),
                      ffiC.gkyl_bc_twistshift_release)
end

function _M.advance(fldDo, fldTar)
   ffiC.gkyl_bc_twistshift_advance(_M._zero, fldDo._zero, fldTar._zero)
end

local p2l = function(valIn, xcIn, dxIn)
   -- Transform a physical coordinate (valIn) to the [-1,1] logical
   -- space in a cell with cell center xcIn and length dxIn.
   return 2.*(valIn-xcIn)/dxIn
end
local dxAxc = function(lims)
   -- Return the length and cell center of the interval [lims.lo, lims.up].
   local lo, up = lims.lo, lims.up
   return up-lo, 0.5*(lo+up)
end

local wrapNum = function(val, lims, pickUpper)
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

local riddersStop = function(eps)   -- Stop Ridder's root-finding if function is < eps.
   return function(x, y, xl, xu, yl, yu)
      if math.abs(y) < eps then return true else return false end
   end
end

function _M.rootFind(funcIn, lims, tol)
   -- Use a Ridder's root finder to find the root of funcIn in the interval
   -- [lims.lo,lims.up] down to a tolerance 'tol'. Return the interval limit
   -- if the function is smaller than the tolerance there. Return nil if the
   -- function does not change sign in the interval (interval doesn't contain the root).
   local funcLo, funcUp = funcIn(lims.lo), funcIn(lims.up)
   if math.abs(funcLo) < tol  then
      return lims.lo
   elseif math.abs(funcUp) < tol  then
      return lims.up
   else
      if funcLo*funcUp < 0. then
         return root.ridders(funcIn, lims.lo, lims.up, riddersStop(tol))
      else
         return nil
      end
   end
end

local yShiftF = function(x) return 1.0 end
function _M.set_yShiftF(shiftFuncIn)
   -- Set the y-shift function that gets called by other functions in this module to
   -- perform the sub-cell integrals. Assuming shiftFuncIn takes (t, xn).
   yShiftF = function(x) return shiftFuncIn(0., {x}) end
end

local domLim = {{lo=0., up=0.}, {lo=0., up=0.}}
function _M.set_domLim(gridIn)
   -- Set table with domain lower and upper boundaries.
   domLim = { {lo=gridIn:lower(1), up=gridIn:upper(1)},
              {lo=gridIn:lower(2), up=gridIn:upper(2)} }
end

local dx = Lin.Vec(2)
function _M.set_dx(gridIn)
   -- Set vector with cell lengths. This gets re-set in the outer loop of preCalcMat
   -- as well (no effect for uniform grids).
   local idx, dxP, dim = {}, {}, gridIn:ndim()
   for d = 1, dim do idx[d] = 1 end
   gridIn:setIndex(idx)
   gridIn:getDx(dxP)
   for d = 1, 2 do dx[d] = dxP[d] end
end

-- Given a y-coordinate of the target cell (yTar), and a y-coordinate
-- of the donor cell (yDo), find the x-coordinate of the point where
-- the yTar-yShift(x) and y=yDo lines intersect. 
--  yTar:    target cell y coordinate.
--  yDo:     donor cell y coordinate.
--  xBounds: search in the interval [xBounds[1],xBounds[2]].
--  yLims:   lower and upper limits of the grid.
-- If y-yShift-yDo=0 has no roots, it is possible that y-yShift intersects a periodic
-- copy of this domain. Check for such cases by looking for the roots of
-- yTar-yShift-(yDo-N*Ly)=0 where Ly is the length of the domain along y and N is an integer.
local findIntersect = function(yTar, yDo, xBounds, yLims)
   local tol = 1.e-13
   local Ly  = yLims.up - yLims.lo
   local function lossF(xIn) return yTar-yShiftF(xIn)-yDo end
   root0 = _M.rootFind(lossF, xBounds, tol)
   if root0 == nil then
      -- Maybe yTar-ySh intersects y=yDo in a periodic copy of the domain. Find roots of
      -- yTar-ySh-(yDo-nP*Ly)=0 where Ly is the y-length of the domain and nP is an integer.
      local nP = {}
      -- Evaluate yTar-ySh at some points in [xBounds.lo,xBounds.up] and obtain potential nP's.
      local stepFrac = 0.1
      local numSteps = math.ceil(1./stepFrac)
      local stepSize = (xBounds.up-xBounds.lo)/numSteps
      for sI = 0, numSteps do
         local xp      = xBounds.lo+sI*stepSize
         local xpShift = yShiftF(xp)
         local exShift = xpShift<0. and ((yTar-xpShift)-yLims.lo)/Ly or (yLims.up-(yTar-xpShift))/Ly
         local newN    = math.sign(xpShift)*math.floor(math.abs(exShift))
         if lume.find(nP,newN)==nil then table.insert(nP,newN) end
      end
      for _, iN in ipairs(nP) do
        local function lossFex(xIn) return yTar-yShiftF(xIn) - (yDo-iN*Ly) end
        root0 = _M.rootFind(lossFex, xBounds, tol)
        if root0 then break end
      end
   end
   return root0
end

local function sign(num)
   -- Sign of a number. Return zero if the number is zero.
   if num > 0 then     return  1
   elseif num < 0 then return -1
   else                return  0
   end
end

local doTarOff = function(xcDoIn, xcTarIn)
   -- y-offset between the donor and the target cell (yDo-yTar), in the direction of the shift.
   --   xcDo:  cell center coordinates of donor cell.
   --   xcTar: cell center coordinates of target cell.
   local shift     = yShiftF(xcDoIn[1])
   local shiftSign = sign(shift)
   local Ly        = domLim[2].up - domLim[2].lo
   -- The idea here is that we keep shifting the donor cell center until it is in a
   -- periodic copy of our domain which overlaps with the shifted target cell center.
   local yDoS, yTarS = xcDoIn[2], xcTarIn[2]-shift
   local keepShifting = true
   while keepShifting do
      local yDoSlo, yDoSup = yDoS-Ly/2., yDoS+Ly/2.
      if yDoSlo <= yTarS and yTarS <= yDoSup then
         keepShifting = false
         break
      else
         yDoS = yDoS - shiftSign*Ly
      end
   end
   return xcTarIn[2] - yDoS
end

-- Given a logical space x coordinate (xi) and a (physical) y-coordinate in the target cell,
-- compute the shifted y-coordinate in the logical space of the donor cell (eta \in [-1,1]).
--   xi:    logical space x coordinate.
--   yTar:  physical y-coordinate in target cell.
--   pmSh:  factor multiplying the y-shift (+/- 1).
--   xcDo:  cell center coordinates of donor cell.
--   xcTar: cell center coordinates of target cell.
--   dx:    cell lengths.
--   pickUpper: boolean indicating if wrapping function should return upper/lower boundary.
local yShiftedLog = function(xi, yTar, pmSh, xcDo, xcTar, pickUpper)
   local xPhys = xcTar[1] + 0.5*dx[1]*xi
   local yS    = yTar-pmSh*yShiftF(xPhys)
   yS = wrapNum(yS, domLim[2], pickUpper)
   local eta   = p2l(yS, xcDo[2], dx[2])
   return eta
end

-- Given an eta coordinate, the donor and target cell centers, compute the
-- corresponding xi. Requires inverting p2l(yTar-yShift(x)).
local yShiftedLogInv = function(ySlog, yTar, xcDoIn, xcTarIn, xLim, pickUpper)
   -- Invert the function p2l(y-yShifted(x)) via root finding. This function is passed to a projection
   -- operation, so that the resulting DG expansion is defined in a subcell of the donor cell.
   local function lossF(xlog)
      return ySlog - yShiftedLog(xlog, yTar, 1., xcDoIn, xcTarIn, pickUpper)
   end
   local tol = 1.e-10
   -- Have to consider extended logical space even though the part of the y-yShift curve
   -- that is not in this cell does not get used in the integral because the fixed
   -- y limit cuts if off. But if we didn't consider extended logical space the shape of the
   -- x-limit function would not come out right.
   local eps = 0.  -- Need wiggle room because it might fail if the root is at the boundary.
   local xlogOut = root.ridders(lossF, xLim.lo-eps, xLim.up+eps, riddersStop(tol))
   return xlogOut
end

-- Permanent data structures needed in nodToModProj1D function and by calls to it.
local projData = {
   xc        = Lin.Vec(1),   -- Cell center coordinate.
   dx        = Lin.Vec(1),   -- Cell length.
   pOrder    = 0,            -- Polynomial order used to project 1D shift function.
   compNodes = {},           -- Node coordinates (computational space).
   numNodes  = 0,
   xPhys_i   = Lin.Vec(1),   -- Physical coordinates at node.
   func_i    = 0,            -- Function evaluated a physical nodes.
   xiLo_eta  = 0,            -- DG coefficients of logical space x-lower limit poly approx.
   xiUp_eta  = 0,            -- DG coefficients of logical space x-upper limit poly approx.
   xiLoL_eta = 0,            -- DG coefficients of logical space x-lower limit poly approx (left).
   xiUpL_eta = 0,            -- DG coefficients of logical space x-upper limit poly approx (left).
   xiLoR_eta = 0,            -- DG coefficients of logical space x-lower limit poly approx (right).
   xiUpR_eta = 0,            -- DG coefficients of logical space x-upper limit poly approx (right).
   etaLo_xi  = 0,            -- DG coefficients of logical space y-lower limit poly approx.
   etaUp_xi  = 0,            -- DG coefficients of logical space y-upper limit poly approx.
   etaLoLo_xi  = 0,            -- DG coefficients of logical space y-lower limit poly approx (lower).
   etaUpLo_xi  = 0,            -- DG coefficients of logical space y-upper limit poly approx (lower).
   etaLoUp_xi  = 0,            -- DG coefficients of logical space y-lower limit poly approx (upper).
   etaUpUp_xi  = 0,            -- DG coefficients of logical space y-upper limit poly approx (upper).
}

function _M.set_projData(shiftPolyOrder)
   -- Create permanent data structures needed by the nodToModProj1D function that
   -- depend on other inputs.
   local compNodes = SerendipityNodes["nodes1xp" .. shiftPolyOrder]
   local numNodes  = #compNodes
   projData.pOrder    = shiftPolyOrder   -- Polynomial order used to project 1D shift function.
   projData.compNodes = compNodes        -- Node coordinates (computational space).
   projData.numNodes  = numNodes
   projData.func_i    = Lin.Mat(numNodes, 1)   -- Function evaluated a physical nodes.
   projData.xiLo_eta  = Lin.Vec(numNodes)      -- DG coefficients of logical space x-lower limit poly approx.
   projData.xiUp_eta  = Lin.Vec(numNodes)      -- DG coefficients of logical space x-upper limit poly approx.
   projData.xiLoL_eta = Lin.Vec(numNodes)      -- DG coefficients of logical space x-lower limit poly approx (left).
   projData.xiUpL_eta = Lin.Vec(numNodes)      -- DG coefficients of logical space x-upper limit poly approx (left).
   projData.xiLoR_eta = Lin.Vec(numNodes)      -- DG coefficients of logical space x-lower limit poly approx (right).
   projData.xiUpR_eta = Lin.Vec(numNodes)      -- DG coefficients of logical space x-upper limit poly approx (right).
   projData.etaLo_xi  = Lin.Vec(numNodes)      -- DG coefficients of logical space y-lower limit poly approx.
   projData.etaUp_xi  = Lin.Vec(numNodes)      -- DG coefficients of logical space y-upper limit poly approx.
   projData.etaLoLo_xi = Lin.Vec(numNodes)      -- DG coefficients of logical space y-lower limit poly approx (lower).
   projData.etaUpLo_xi = Lin.Vec(numNodes)      -- DG coefficients of logical space y-upper limit poly approx (lower).
   projData.etaLoUp_xi = Lin.Vec(numNodes)      -- DG coefficients of logical space y-lower limit poly approx (upper).
   projData.etaUpUp_xi = Lin.Vec(numNodes)      -- DG coefficients of logical space y-upper limit poly approx (upper).
end

-- Functions used by nodToModProj1D (copied from EvalOnNodes).
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
ffi.cdef([[ void nodToMod(double* fN, int numNodes, int numVal, int ndim, int p, double* fM); ]])
local evalFunc   = loadstring(evalFuncTempl { M = 1 } )()
local compToPhys = loadstring(compToPhysTempl {NDIM = 1} )()

local nodToModProj1D = function(funcIn, limIn, fldOut)
   -- Project 'funcIn' onto 1D DG basis in the interval [limIn.lo, limIn.up]
   -- evaluating at nodes and transforming to modal representation.
   --   funcIn: 1D scalar function to be projected.
   --   limIn:  limits of the interval in which to project the function.
   --   fldOut: DG field output.

   projData.dx[1], projData.xc[1] = dxAxc(limIn)
   for i = 1, projData.numNodes do
      compToPhys(projData.compNodes[i], projData.dx, projData.xc, projData.xPhys_i)
      evalFunc(0.0, projData.xPhys_i, funcIn, projData.func_i[i])
   end
   ffi.C.nodToMod(projData.func_i:data(), projData.numNodes, 1, 1, projData.pOrder, fldOut:data())
end

local subCellInt_xLimDG = function(xLimFuncs, yBounds, xIdxIn, offDoTar, yShPtrIn, jDo)
   -- Perform a sub-cell integral with xi-limits that are functions of eta given by a
   -- DG polynomial representation. Here xi-eta mean the logical coordinates of the donor cell.
   -- The functions xLimFuncs define these xi-limits in the eta-space of the
   -- donor cell. Then we project these functions onto a 1D DG basis. Since the integral
   -- may not span the whole eta-space of the donor cell (dySub<2), the 1D
   -- projection is done treating the partial eta-space of the donor cell as a
   -- 'physical' space on which the xi-limit polynomial is defined, which has a different
   -- relationship to its own logical coordinate than that between the logical and physical
   -- coordinate of the donor cell. This means that the difference between the logical
   -- coordinates of the xi-limits and the donor field needs to be accounted for. We do that
   -- using dySub and ycSub in the kernel.
   --   xLimFuncs: functions defining lower/upper x-limits of integral.
   --   yBounds:   lower/upper logical y-limits of the integral.
   --   xIdx:      current x-index.
   --   offDoTar:  physical y-offset between donor and target cells.
   --   yShPtr:    pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local dy       = dx[2]
   local dySub, _ = dxAxc(yBounds)
   if dySub > 0 then
      nodToModProj1D(xLimFuncs.lo, yBounds, projData.xiLo_eta)
      nodToModProj1D(xLimFuncs.up, yBounds, projData.xiUp_eta)
      ffiC.gkyl_bc_twistshift_integral_xlimdg(_M._zero, 1., projData.xiLo_eta:data(), projData.xiUp_eta:data(), yBounds.lo, yBounds.up, dy, offDoTar, yShPtrIn, xIdxIn[1], jDo-1) -- need the do index

   end
end

local subCellInt_xLoLimDG = function(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, pickUpper, yShPtrIn, jDo)
   -- Perform subcell integral with variable lower x-limit.
   --   yTar:      yTar-yShift gives the curve defining the lower x-limit.
   --   xiBounds:  logical x range in which to invert yTar-yShift(x).
   --   etaBounds: logical y limits of the integral.
   --   xIdxIn:    current x-index.
   --   xcDo:      cell center coordinates of donor cell.
   --   xcTar:     cell center coordinates of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   --   yShPtrIn:  pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local ycOff   = doTarOff(xcDoIn, xcTarIn)
   local xiLims  = {lo = function(t,xn)
                       return yShiftedLogInv(xn[1], yTar, xcDoIn, xcTarIn, xiBounds, pickUpper)
                    end,
                    up = function(t,xn) return 1.0 end}
   subCellInt_xLimDG(xiLims, etaBounds, xIdxIn, ycOff, yShPtrIn, jDo)
end

local subCellInt_xUpLimDG = function(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, pickUpper, yShPtrIn, jDo)
   -- Perform subcell integral with variable upper x-limit.
   --   yTar:      yTar-yShift gives the curve defining the upper x-limit.
   --   xiBounds:  xi range in which to invert yTar-yShift(x).
   --   etaBounds: eta y limits of the integral.
   --   xIdxIn:    current x-index.
   --   xcDo:      cell center coordinates of donor cell.
   --   xcTar:     cell center coordinates of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   --   yShPtrIn:  pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local ycOff    = doTarOff(xcDoIn, xcTarIn)
   local xiLims = {lo = function(t,xn) return -1.0 end,
                   up = function(t,xn)
                      return yShiftedLogInv(xn[1], yTar, xcDoIn, xcTarIn, xiBounds, pickUpper)
                   end}
   subCellInt_xLimDG(xiLims, etaBounds, xIdxIn, ycOff, yShPtrIn, jDo)
end

local subCellInt_yLimDG = function(xBounds, yLimFuncs, xIdxIn, offDoTar, yShPtrIn, jDo)
   -- Perform a sub-cell integral with eta-limits that are functions of xi given by a
   -- DG polynomial representation. Here xi-eta mean the logical coordinates of the donor cell.
   -- The functions yLimFuncs define these eta-limits in the xi-space of the
   -- donor cell. Then we project these functions onto a 1D DG basis. Since the integral
   -- may not span the whole xi space of the donor cell (dxSub<2), the 1D
   -- projection is done treating the partial xi space of the donor cell as a
   -- 'physical' space on which the eta-limit polynomial is defined, which has a different
   -- relationship to its own logical coordinate than that between the logical and physical
   -- coordinate of the donor cell. This means that the difference between the logical
   -- coordinates of the eta-limits and the donor field needs to be accounted for. We do that
   -- using dxSub and xcSub in the kernel.
   --   xBounds:   lower/upper logical x-limits of the integral.
   --   yLimFuncs: functions defining lower/upper y-limits of integral.
   --   xIdx:      current x-index.
   --   offDoTar:  physical y-offset between donor and target cells.
   --   yShPtr:    pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local dy       = dx[2]
   local dxSub, _ = dxAxc(xBounds)
   if dxSub > 0 then
      nodToModProj1D(yLimFuncs.lo, xBounds, projData.etaLo_xi)
      nodToModProj1D(yLimFuncs.up, xBounds, projData.etaUp_xi)
      ffiC.gkyl_bc_twistshift_integral_ylimdg(_M._zero, 1., xBounds.lo, xBounds.up, projData.etaLo_xi:data(), projData.etaUp_xi:data(), dy, offDoTar, yShPtrIn, xIdxIn[1], jDo-1)


   end
end

local subCellInt_yLoLimDG = function(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, pickUpper, yShPtrIn, jDo)
   -- Perform subcell integral with variable lower eta-limit.
   --   yTar:      yTar-yShift gives the curve defining the lower x-limit.
   --   xiBounds:  logical x range in which to invert yTar-yShift(x).
   --   etaBounds: logical y limits of the integral.
   --   xIdxIn:    current x-index.
   --   xcDo:      cell center coordinates of donor cell.
   --   xcTar:     cell center coordinates of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   --   yShPtrIn:  pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local ycOff   = doTarOff(xcDoIn, xcTarIn)
   local etaLims = {lo = function(t,xn)
                        return yShiftedLog(xn[1], yTar, 1, xcDoIn, xcTarIn, pickUpper)
                    end,
                    up = function(t,xn) return 1.0 end}
   subCellInt_yLimDG(xiBounds, etaLims, xIdxIn, ycOff, yShPtrIn, jDo)
end

local subCellInt_yUpLimDG = function(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, pickUpper, yShPtrIn, jDo)
   -- Perform subcell integral with variable upper y-limit.
   --   yTar:      yTar-yShift gives the curve defining the upper x-limit.
   --   xiBounds:  xi range in which to invert yTar-yShift(x).
   --   etaBounds: eta y limits of the integral.
   --   xIdxIn:    current x-index.
   --   xcDo:      cell center coordinates of donor cell.
   --   xcTar:     cell center coordinates of target cell.
   --   pickUpper: boolean to indicate if wrapping function should return upper instead of lower boundary.
   --   yShPtrIn:  pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local ycOff   = doTarOff(xcDoIn, xcTarIn)
   local etaLims = {lo = function(t,xn) return -1.0 end,
                    up = function(t,xn)
                       return yShiftedLog(xn[1], yTar, 1, xcDoIn, xcTarIn, pickUpper)
                    end}
   subCellInt_yLimDG(xiBounds, etaLims, xIdxIn, ycOff, yShPtrIn, jDo)
end

local subCellInt_trapezoid = function(xiBounds, etaBounds, yTar, xIdxIn, xcDoIn, xcTarIn, pickUpper, yShPtrIn, jDo)
   -- Perform integral over a trapezoidal subcell region.
   --   xiBounds:  2x2 table with the bounds of the xi lower and upper limits.
   --   etaBounds: 2 element table with the logical y-bounds of the logical x lower and upper limits.
   --   yTar:      y_{j_tar-/+1/2}-yShift defines the lower/upper x limits.
   --   xIdx:      current x-index.
   --   xcDo:      cell center coordinates of donor cell.
   --   xcTar:     cell center coordinates of target cell.
   --   pickUpper: boolean indicating if wrapping function should return upper/lower boundary.
   --   yShPtrIn:  pointers to y-shift field data in current cell.
   --   jDo: y index of donor cell
   local ycOff = doTarOff(xcDoIn, xcTarIn)
   -- Functions describing the lower and upper limits in logical [-1,1] x-space.
   local xiLims = {
      lo = function(t,xn)
         return yShiftedLogInv(xn[1], yTar.lo, xcDoIn, xcTarIn, xiBounds.lo, pickUpper)
      end,                                                        
      up = function(t,xn)                                         
         return yShiftedLogInv(xn[1], yTar.up, xcDoIn, xcTarIn, xiBounds.up, pickUpper)
      end
   }
   subCellInt_xLimDG(xiLims, etaBounds, xIdxIn, ycOff, yShPtrIn, jDo)
end

local subCellInt_sN = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sN subcell integral: the subcell region is a 4-sided
   -- trapezoid. Both y_{j_tar-/+1/2}-yShift intersect y_{j_do-/+1/2}.
   --   x_pq:         intersections of y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2}.
   --   xIdxIn:       current x-index.
   --   xcDo:         cell center of donor cell.
   --   xcTar:        cell center of target cell.
   --   limTar:       cell limits of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtrIn:     pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local etaBounds = {lo=-1., up=1.}   -- Logical-space y-limits of the subcell integral.

   printDebug("    scenario sN")
   local yTar      = {lo=nil, up=nil}  -- yTar-yShift defines the x-limits of the subcell integral.
   -- Lo/up limits of x-logical space of lo/up limits of integral.
   local xiBounds = {lo={lo=nil, up=nil}, up={lo=nil, up=nil}}
   -- Establish if the y_{j_tar-/+1/2}-yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=y_{j_do-1/2} crosses y_{j_tar-/+1/2}-yShift.
   if x_pq.loTar.loDo < x_pq.upTar.loDo then   -- Monotonically increasing yShift. Trapezoid slanted left.
      yTar.lo, yTar.up = limTar[2].lo, limTar[2].up
      xiBounds.lo.lo, xiBounds.lo.up = p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])
      xiBounds.up.lo, xiBounds.up.up = p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])
   else                              -- Monotonically decreasing yShift. Trapezoid slanted right.
      yTar.lo, yTar.up = limTar[2].up, limTar[2].lo
      xiBounds.lo.lo, xiBounds.lo.up = p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])
      xiBounds.up.lo, xiBounds.up.up = p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])
   end

   subCellInt_trapezoid(xiBounds, etaBounds, yTar, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_siORsii = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario si or sii subcell integral: the subcell region is 3-sided, abutting
   -- one of the upper corners of the donor cell. Only y_{j_tar-1/2}-yShift intersects
   -- y_{j_do+1/2}.
   --   x_pq:         points where y_{j_tar-/+1/2}-yShift intersect y_{j_do-/+1/2}.
   --   xIdxIn:       current x-index.
   --   xcDo:         cell center of donor cell.
   --   xcTar:        cell center of target cell.
   --   limTar:       cell limits of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtrIn:     pointer to y-shift field in current cell.
   --   jDo: y index of donor cell

   -- Two options below:
--   -- Option A: integral with variable x-limits.
--   local yShiftUp = yShiftF(limTar[1].up)
--   local yTar     = limTar[2].lo
--   if -yShiftUp >= -yShiftF(x_pq.loTar.upDo) then  -- Check if yShift increases or decreases.
--      printDebug("    scenario si")
--      -- Scenario si.
--      local yLogBounds = {lo=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
--      local xLogBounds = {lo=-1, up=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])}
--      subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, tsMatVecsIn)
--   else
--      -- Scenario sii.
--      printDebug("    scenario sii")
--      local yLogBounds = {lo=p2l(wrapNum(limTar[2].lo-yShiftUp,domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
--      local xLogBounds = {lo=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), up=1.}
--      subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, tsMatVecsIn)
--   end

   -- Option B: integral with variable y-limits.
   local yShiftUp = yShiftF(limTar[1].up)
   local yTar     = limTar[2].lo
   local etaBounds, xiBounds
   if -yShiftUp >= -yShiftF(x_pq.loTar.upDo) then  -- Check if yShift increases or decreases.
      printDebug("    scenario si")
      -- Scenario si.
      etaBounds = {lo=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
      xiBounds  = {lo=-1, up=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])}
   else
      -- Scenario sii.
      printDebug("    scenario sii")
      etaBounds = {lo=p2l(wrapNum(limTar[2].lo-yShiftUp,domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
      xiBounds  = {lo=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), up=1.}
   end
   subCellInt_yLoLimDG(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_siiiORsiv = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario siii or siv subcell integral: the subcell region is 3-sided, abutting
   -- one of the lower corners of the donor cell. Only y_{j + 1/2}-yShift intersects
   -- y_{j+m - 1/2} (the upper y-boundary of donor cell).
   --   x_pq:         intersections of y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2}.
   --   xIdxIn:       current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtrIn:     pointer to y-shift field in current cell.
   --   jDo: y index of donor cell

   -- Two options below:
--   -- Option A: integral with variable x-limits.
--   local yShiftUp = yShiftF(limTar[1].up)
--   local yTar     = limTar[2].up
--   if -yShiftUp <= -yShiftF(x_pq.upTar.loDo) then
--      -- Scenario siii.
--      printDebug("    scenario siii")
--      local yLogBounds = {lo=-1., up=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
--      local xLogBounds = {lo=-1., up=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])}
--      subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, tsMatVecsIn)
--   else
--      printDebug("    scenario siv")
--      -- Scenario siv.
--      local yLogBounds = {lo=-1., up=p2l(wrapNum(limTar[2].up-yShiftUp,domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
--      local xLogBounds = {lo=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), up=1.}
--      subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, tsMatVecsIn)
--   end

   -- Option B: integral with variable y-limits.
   local yShiftUp = yShiftF(limTar[1].up)
   local yTar     = limTar[2].up
   local etaBounds, xiBounds
   if -yShiftUp <= -yShiftF(x_pq.upTar.loDo) then
      -- Scenario siii.
      printDebug("    scenario siii")
      etaBounds = {lo=-1., up=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
      xiBounds  = {lo=-1., up=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])}
   else
      printDebug("    scenario siv")
      -- Scenario siv.
      etaBounds = {lo=-1., up=p2l(wrapNum(limTar[2].up-yShiftUp,domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
      xiBounds  = {lo=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), up=1.}
   end
   subCellInt_yUpLimDG(yTar, xiBounds, etaBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_six = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario six subcell integral.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local xLogBounds = {lo=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])}
   local yLogBounds = {lo=-1., up=1.}
   local yTar       = limTar[2].lo    -- y_{j_tar-/+1/2} to shift in y.
   subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_sx = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sx subcell integral.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local xLogBounds = {lo=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])}
   local yLogBounds = {lo=-1., up=1.}
   local yTar       = limTar[2].lo    -- y_{j_tar-/+1/2} to shift in y.
   subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_sxi = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sxi subcell integral.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+ 1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local xLogBounds = {lo=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])}
   local yLogBounds = {lo=-1., up=1.}
   local yTar       = limTar[2].up    -- y_{j_tar-/+1/2} to shift in y.
   subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_sxii = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sxii subcell integral.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local xLogBounds = {lo=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])}
   local yLogBounds = {lo=-1., up=1.}
   local yTar       = limTar[2].up
   subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
end

local subCellInt_svORsvi = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sv or svi subcell integral: the subcell region is 5-sided, abutting
   -- one of the lower corners of the donor cell. Only y_{j_tar+1/2}-yShift does not intersect
   -- y_{j_do-1/2}.
   --   x_pq:         intersections of y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   -- Establish if the y_{j_tar-/+1/2}-yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=y_{j_do+1/2} crosses y_{j_tar-/+1/2}-yShift.
   if x_pq.loTar.upDo > x_pq.upTar.upDo then   -- Scenario sv.
      local xLogLoLo, xLogLoUp = -1., p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])
      printDebug("    scenario sv")
      if (xLogLoUp-xLogLoLo) > 1.e-13 then
         -- Add scenario sN-like contribution.
         local yTars      = {lo=limTar[2].up, up=limTar[2].lo}
         local xLogBounds = {lo={lo=xLogLoLo, up=xLogLoUp},
                             up={lo=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])}}
         local yLogBounds = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
         subCellInt_trapezoid(xLogBounds, yLogBounds, yTars, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
         -- Add scenario siii-like contribution.
         local yTar    = limTar[2].lo
         xLogBounds    = {lo=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1])}
         yLogBounds.up = yLogBounds.lo    -- Keep the order of yLogBounds here.
         yLogBounds.lo = -1.
         subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
      else   -- This is better handled as a scenario six.
         printDebug("      as six")
         subCellInt_six(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
      end
   else   -- Scenario svi.
      local xLogUpLo, xLogUpUp = p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), 1.
      printDebug("    scenario svi")
      if (xLogUpUp-xLogUpLo) > 1.e-13 then
         -- Add scenario sN-like contribution.
         local yTars      = {lo=limTar[2].lo, up=limTar[2].up}
         local xLogBounds = {lo={lo=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])},
                             up={lo=xLogUpLo, up=xLogUpUp}}
         local yLogBounds = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].up),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
         subCellInt_trapezoid(xLogBounds, yLogBounds, yTars, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
         -- Add scenario siv-like contribution.
         local yTar    = limTar[2].lo
         xLogBounds    = {lo=p2l(x_pq.loTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])}
         yLogBounds.up = yLogBounds.lo    -- Keep the order of yLogBounds here.
         yLogBounds.lo = -1.
         subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
      else    -- This is better handled as a scenario sx.
         printDebug("      as sx")
         subCellInt_sx(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
      end
   end
end

local subCellInt_sviiORsviii = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario svii or sviii subcell integral.
   -- The subcell region is 5-sided, abutting one of the upper corners of the donor cell.
   -- Only y_{j_tar-1/2}-yShift does not intersects y_{j_do+1/2}.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+1/2}.
   --   xIdx:         current x-index.
   --   xcDo:         cell center coordinates of donor cell.
   --   xcTar:        cell center coordinates of target cell.
   --   limTar:       cell boundaries along each direction of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtr:       pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   -- Establish if the y_{j_tar-/+1/2}-yShift defines the lower or upper limit of the
   -- x-integral by comparing where y=y_{j_do-1/2} crosses y_{j_tar-/+1/2}-yShift.
   if x_pq.loTar.loDo < x_pq.upTar.loDo then   -- Scenario svii.
      local xLimLoLo, xLimLoUp = -1., p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])
      printDebug("    scenario svii")
      if (xLimLoUp-xLimLoLo) > 1.e-13 then
         -- Add scenario sN-like contribution.
         local yTar       = {lo=limTar[2].lo, up=limTar[2].up}
         local xLogBounds = {lo={lo=xLimLoLo, up=xLimLoUp},
                             up={lo=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])}}
         local yLogBounds = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
         subCellInt_trapezoid(xLogBounds, yLogBounds, yTar, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
         -- Add scenario si-like contribution.
         local yTar    = limTar[2].up
         xLogBounds    = {lo=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1])}
         yLogBounds.lo = yLogBounds.up    -- Keep the order of yLogBounds here.
         yLogBounds.up = 1.
         subCellInt_xUpLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
      else   -- This is better handled as a scenario sxi.
         printDebug("      as sxi")
         subCellInt_sxi(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
      end
   else   -- Scenario sviii.
      local xLimUpLo, xLimUpUp = p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), 1.
      printDebug("    scenario sviii")
      if (xLimUpUp-xLimUpLo) > 1.e-13 then
         -- Add scenario sN-like contribution.
         local yTar       = {lo=limTar[2].up, up=limTar[2].lo}
         local xLogBounds = {lo={lo=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])},
                             up={lo=xLimUpLo, up=xLimUpUp}}
         local yLogBounds = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].up),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
         subCellInt_trapezoid(xLogBounds, yLogBounds, yTar, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
         -- Add scenario sii-like contribution.
         local yTar       = limTar[2].up
         local xLogBounds = {lo=p2l(x_pq.upTar.loDo,xcDoIn[1],dx[1]), up=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])}
         yLogBounds.lo    = yLogBounds.up    -- Keep the order of yLogBounds here.
         yLogBounds.up    = 1.
         subCellInt_xLoLimDG(yTar, xLogBounds, yLogBounds, xIdxIn, xcDoIn, xcTarIn, atUpperYcell, yShPtrIn, jDo)
      else   -- This is better handled as a scenario sxii.
         printDebug("      as sxii")
         subCellInt_sxii(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
      end
   end
end

local subCellInt_sxiiiORsxiv = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario svxiii or sxiv subcell integral.
   -- The sub-cell region consists of the integral over the full cell minus the integral
   -- over two opposing corners.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift of y_{j_do-/+1/2} (upper/lower y-boundaries of donor cell).
   --   xIdxIn:       current x-index.
   --   xcDo:         cell center of donor cell.
   --   xcTar:        cell center of target cell.
   --   limTar:       cell limits of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtrIn:     pointer to y-shift field in current cell.
   --   jDo: y index of donor cell

   local ycOff = doTarOff(xcDoIn, xcTarIn)
   -- Add the contribution from integrating over the whole cell.
   ffiC.gkyl_bc_twistshift_integral_fullcelllimdg(_M._zero, dx[2], ycOff, yShPtrIn, xIdxIn[1], jDo-1)


   -- Two options below:
--   -- Option A: integral with variable x-limits.
--   local yTarL, yTarR
--   local xLogBoundsL, xLogBoundsR = {}, {}
--   local yLogBoundsL, yLogBoundsR = {}, {}
--   if -yShiftF(limTar[1].lo) < -yShiftF(limTar[1].up) then   -- Scenario sxiii.
--      printDebug("    scenario sxiii")
--      -- Scenario si-like part.
--      yTarL       = limTar[2].up
--      xLogBoundsL = {lo=-1., up=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])}
--      yLogBoundsL = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
--      -- Scenario siv-like part.
--      yTarR       = limTar[2].lo
--      xLogBoundsR = {lo=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), up=1.}
--      yLogBoundsR = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].up),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
--   else   -- Scenario sxiv.
--      printDebug("    scenario sxiv")
--      -- Scenario sii-like part.
--      yTarR       = limTar[2].up
--      xLogBoundsR = {lo=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), up=1.}
--      yLogBoundsR = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].up),domLim[2],atUpperYcell), xcDoIn[2], dx[2]), up=1.}
--      -- Scenario siii-like part.
--      yTarL       = limTar[2].lo
--      xLogBoundsL = {lo=-1., up=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])}
--      yLogBoundsL = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
--   end
--   -- Project x-limits of left subcell region to subtract.
--   local xLogLimsL = {lo = function(t,xn) return -1.0 end,
--                      up = function(t,xn)
--                         return yShiftedLogInv(xn[1], yTarL, xcDoIn, xcTarIn, xLogBoundsL, atUpperYcell)
--                      end}
--   local dySubL, _ = dxAxc(yLogBoundsL)
--   if (xLogBoundsL.up-xLogBoundsL.lo) > 1.e-13 and dySubL > 1.e-13 then
--      nodToModProj1D(xLogLimsL.lo, yLogBoundsL, projData.xiLoL_eta)
--      nodToModProj1D(xLogLimsL.up, yLogBoundsL, projData.xiUpL_eta)
--      -- Subtract the left corner subcell region.
--      intSubXlimDG(-1., projData.xiLoL_eta:data(), projData.xiUpL_eta:data(), yLogBoundsL.lo, yLogBoundsL.up,
--                   dx[2], ycOff, yShPtrIn, xIdxIn[1], 0)
--   end
--   -- Project x-limits of right subcell region to subtract.
--   local xLogLimsR = {lo = function(t,xn)
--                         return yShiftedLogInv(xn[1], yTarR, xcDoIn, xcTarIn, xLogBoundsR, atUpperYcell)
--                      end,
--                      up = function(t,xn) return 1.0 end}
--   local dySubR, _ = dxAxc(yLogBoundsR)
--   if (xLogBoundsR.up-xLogBoundsR.lo) > 1.e-13 and dySubR > 1.e-13 then
--      nodToModProj1D(xLogLimsR.lo, yLogBoundsR, projData.xiLoR_eta)
--      nodToModProj1D(xLogLimsR.up, yLogBoundsR, projData.xiUpR_eta)
--      -- Subtract the right corner subcell region.
--      intSubXlimDG(-1., projData.xiLoR_eta:data(), projData.xiUpR_eta:data(), yLogBoundsR.lo, yLogBoundsR.up,
--                   dx[2], ycOff, yShPtrIn, xIdxIn[1], 0)
--   end


   -- Option B: integral with variable y-limits.
   local yTarLo, yTarUp
   local xiBoundsLo, xiBoundsUp  = {}, {}
   local etaBoundsLo, etaBoundUp = {}, {}
   if -yShiftF(limTar[1].lo) < -yShiftF(limTar[1].up) then   -- Scenario sxiii.
      printDebug("    scenario sxiii")
      -- Scenario si-like part.
      yTarUp      = limTar[2].up
      xiBoundsUp  = {lo=-1., up=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1])}
      etaBoundsUp = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2]), up=1.}
      -- Scenario siv-like part.
      yTarLo      = limTar[2].lo
      xiBoundsLo  = {lo=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1]), up=1.}
      etaBoundsLo = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].up),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
   else   -- Scenario sxiv.
      printDebug("    scenario sxiv")
      -- Scenario sii-like part.
      yTarUp      = limTar[2].up
      xiBoundsUp  = {lo=p2l(x_pq.upTar.upDo,xcDoIn[1],dx[1]), up=1.}
      etaBoundsUp = {lo=p2l(wrapNum(limTar[2].up-yShiftF(limTar[1].up),domLim[2],atUpperYcell), xcDoIn[2], dx[2]), up=1.}
      -- Scenario siii-like part.
      yTarLo      = limTar[2].lo
      xiBoundsLo  = {lo=-1., up=p2l(x_pq.loTar.loDo,xcDoIn[1],dx[1])}
      etaBoundsLo = {lo=-1., up=p2l(wrapNum(limTar[2].lo-yShiftF(limTar[1].lo),domLim[2],atUpperYcell),xcDoIn[2],dx[2])}
   end

   local etaLimsUp = {lo = function(t,xn)
                         return yShiftedLog(xn[1], yTarUp, 1, xcDoIn, xcTarIn, atUpperYcell)
                      end,
                      up = function(t,xn) return 1.0 end}
   local dxSubUp, _ = dxAxc(xiBoundsUp)
   if (etaBoundsUp.up-etaBoundsUp.lo) > 1.e-13 and dxSubUp > 1.e-13 then
      nodToModProj1D(etaLimsUp.lo, xiBoundsUp, projData.etaLoUp_xi)
      nodToModProj1D(etaLimsUp.up, xiBoundsUp, projData.etaUpUp_xi)
      ffiC.gkyl_bc_twistshift_integral_ylimdg(_M._zero,-1., xiBoundsUp.lo, xiBoundsUp.up, projData.etaLoUp_xi:data(), projData.etaUpUp_xi:data(),
                   dx[2], ycOff, yShPtrIn, xIdxIn[1], jDo-1)

   end

   local etaLimsLo = {lo = function(t,xn) return -1.0 end,
                      up = function(t,xn)
                         return yShiftedLog(xn[1], yTarLo, 1, xcDoIn, xcTarIn, atUpperYcell)
                      end}
   local dxSubLo, _ = dxAxc(xiBoundsLo)
   if (etaBoundsLo.up-etaBoundsLo.lo) > 1.e-13 and dxSubLo > 1.e-13 then
      nodToModProj1D(etaLimsLo.lo, xiBoundsLo, projData.etaLoLo_xi)
      nodToModProj1D(etaLimsLo.up, xiBoundsLo, projData.etaUpLo_xi)
      ffiC.gkyl_bc_twistshift_integral_ylimdg(_M._zero, -1., xiBoundsLo.lo, xiBoundsLo.up, projData.etaLoLo_xi:data(), projData.etaUpLo_xi:data(), dx[2], ycOff, yShPtrIn, xIdxIn[1], jDo-1)

   end
end

local subCellInt_sxvORsxvi = function(x_pq, xIdxIn, xcDoIn, xcTarIn, limDo, limTar, atUpperYcell, yShPtrIn, jDo)
   -- Perform scenario sxv or sxvi subcell integral, using fixed x-limits and variable y limits.
   --   x_pq:         intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
   --   xIdxIn:       current x-index.
   --   xcDo:         cell center of donor cell.
   --   xcTar:        cell center of target cell.
   --   limTar:       cell limits of target cell.
   --   atUpperYcell: boolean to indicate if this donor is the upper Y cell.
   --   yShPtrIn:     pointer to y-shift field in current cell.
   --   jDo: y index of donor cell
   local xiBounds = {lo=-1., up=1.}  -- Limits of xi integral.
   local etaLims  = {}               -- Table of functions definting the limits of eta integral.
   local yShiftLoXc = wrapNum(limTar[2].lo-yShiftF(xcDoIn[1]),domLim[2],atUpperYcell)
   if (limDo[2].lo <= yShiftLoXc) and (yShiftLoXc <= limDo[2].up) then   -- Scenario sxv.
      printDebug("    scenario sxv")
      local yTar = limTar[2].lo
      etaLims = {lo = function(t,xn)
                         return yShiftedLog(xn[1], yTar, 1, xcDoIn, xcTarIn, atUpperYcell)
                      end,
                 up = function(t,xn) return 1.0 end}
   else   -- Scenario sxvi.
      printDebug("    scenario sxvi")
      local yTar = limTar[2].up
      etaLims = {lo = function(t,xn) return -1.0 end,
                 up = function(t,xn)
                         return yShiftedLog(xn[1], yTar, 1, xcDoIn, xcTarIn, atUpperYcell)
                      end}
   end
   nodToModProj1D(etaLims.lo, xiBounds, projData.etaLo_xi)
   nodToModProj1D(etaLims.up, xiBounds, projData.etaUp_xi)
   local ycOff = doTarOff(xcDoIn, xcTarIn)  -- Offset between cell centers along y.
   ffiC.gkyl_bc_twistshift_integral_ylimdg(_M._zero, 1., xiBounds.lo, xiBounds.up, projData.etaLo_xi:data(), projData.etaUp_xi:data(),
                dx[2], ycOff, yShPtrIn, xIdxIn[1], jDo-1)

end


function _M.getDonors(grid, yShift, yShBasis)
   -- Identify the donor cells for each target cell. This only considers
   -- a single 2D plane and returns 2D indices.
   --   grid:     grid object on which donor/target fields is defined.
   --   yShift:   1D field with the y-shift.
   --   yShBasis: 1D basis of yShift.
   
   local doCells = {}   -- Donor cells, one table for each target cell.
   for i = 1, grid:numCells(1) do
      doCells[i] = {}
      for j = 1, grid:numCells(2) do doCells[i][j] = {} end
   end

   -- Determine donor cells by taking points just inside the boundary of the target cell,
   -- shifting them, and finding out which cells own the shifted points. 

   local stepFrac  = {0.1, 0.1}  -- Step taken around boundary, as fraction of cell length.
   local deltaFrac = 1.e-9       -- Distance away from the boundary, as fraction of cell length.
   local numSteps  = {}
   for d = 1,2 do numSteps[d] = math.ceil(1./stepFrac[d]) end

   local yShIndexer, yShItr = yShift:genIndexer(), yShift:get(1)

   local yShNumB    = yShift:numComponents()
   local yShBasisEv = Lin.Vec(yShNumB)

   local cellLim, stepSz, delta  = {{lo=0., up=0.},{lo=0., up=0.}}, {0., 0.}, {0., 0.} 
   local evPoint                 = {0., 0.}
   local newP, newIdx            = {0., 0.}, {0, 0}
   local idxShifted, searchPoint, idxP = {}, {}, {}
   for d = 1, grid:ndim() do
      idxShifted[d], idxP[d] = 1, 1
      searchPoint[d] = grid:lower(d)+deltaFrac*grid:dx(d)
   end

   local xyRange = grid:localRange():selectFirst(2)
   for idx in xyRange:rowMajorIter() do
      grid:setIndex(idx)
      
      local doCellsC = doCells[idx[1]][idx[2]]  -- Donor cells of current target cell.

      for d = 1,2 do
         cellLim[d].lo, cellLim[d].up = grid:cellLowerInDir(d), grid:cellUpperInDir(d)
         stepSz[d] = grid:dx(d)/numSteps[d]
         delta[d]  = deltaFrac*grid:dx(d)
      end

      yShift:fill(yShIndexer(idx), yShItr)

      for dC = 1,2 do            -- dC=1: y=const, dC=2: x=const (boundaries).
         for xS = 1,2 do         -- xS=1: lower, xS=2 upper (boundary).

            -- Search first shifted point. Use pickLower=false in findCell unless
            -- searching for points along a x=const line near the upper y-boundary.
            for d = 1,2 do evPoint[d] = cellLim[d].lo+delta[d] end
            evPoint[dC] = evPoint[dC]+(xS-1)*(grid:dx(dC)-2.*delta[dC])

            -- Evaluate yShift at this point.
            local p2l_eval= Lin.Vec(1)
            p2l_eval[1] = p2l(evPoint[1],grid:cellCenterInDir(1),grid:dx(1))
            yShBasis:evalBasis(p2l_eval:data(), yShBasisEv:data())
            local yShEv = 0.
            for k = 1,yShNumB do yShEv = yShEv + yShItr[k]*yShBasisEv[k] end

            -- Find the index of the cell that owns the shifted point.
            for d = 1,2 do idxShifted[d] = nil end
            searchPoint[1], searchPoint[2] = evPoint[1], wrapNum(evPoint[2]-yShEv,domLim[2],false)
            local chooseLo = (dC==2 and xS==2) and true or false
            idxP[1], idxP[2] = idx[1], nil
            grid:findCell(searchPoint,idxShifted,chooseLo,idxP)
            newIdx = {idxShifted[1], idxShifted[2]}
            if lume.findTable(doCellsC,newIdx)==nil then table.insert(doCellsC,newIdx) end

            -- Search other shifted points along this line.
            local stepDim = dC % 2+1
            newP     = {nil,nil}
            newP[dC] = evPoint[dC]
            for sI = 1, numSteps[stepDim] do

               newP[stepDim] = evPoint[stepDim] + sI*stepSz[stepDim]
               if sI==numSteps[stepDim] then newP[stepDim]=newP[stepDim]-2.*delta[stepDim] end

               -- Evaluate yShift at this point.
               p2l_eval[1] = p2l(newP[1],grid:cellCenterInDir(1),grid:dx(1))
               yShBasis:evalBasis(p2l_eval:data(), yShBasisEv:data())
               yShEv = 0.
               for k = 1,yShNumB do yShEv = yShEv + yShItr[k]*yShBasisEv[k] end

               -- Find the index of the cell that owns the shifted point.
               for d = 1,2 do idxShifted[d] = nil end
               searchPoint[1], searchPoint[2] = newP[1], wrapNum(newP[2]-yShEv,domLim[2],false)
               idxP[1], idxP[2] = idx[1], nil
               grid:findCell(searchPoint,idxShifted,true,idxP)
               newIdx = {idxShifted[1], idxShifted[2]}
               if lume.findTable(doCellsC,newIdx)==nil then table.insert(doCellsC,newIdx) end
            end

         end
      end

   end

   return doCells
end

-- Count nils in the table of intersection points and provide
-- their indices. Also give the indices of the non-nil elements.
local missingInterPt = function(x_pq)
  local count, nilIdxs, nonNilIdxs = 0, {}, {}
  for i, _ in pairs(x_pq) do
     for _, j in ipairs({"loDo","upDo"}) do
        if x_pq[i][j] == nil then
           count = count+1
           table.insert(nilIdxs,{i,j})
        else
           table.insert(nonNilIdxs,{i,j})
        end
     end
  end
  return count, nilIdxs, nonNilIdxs
end

function _M.preCalcMat(grid, yShift, doCells)
   -- Pre-calculate the matrices that multiply DG coefficients of each donor cell to obtain their
   -- contribution to each target cell. Based on weak equality between donor and target fields.
   --   grid:      grid object on which donor/target fields is defined.
   --   yShift:    1D field with the y-shift.
   --   doCells:   donor cell 2D indices.

   local dim = grid:ndim()

   local localRange = grid:localRange()

   local yIdxTar = 1   -- For positive(negative) yShift idx=1(last) might be better, but ideally it shouldn't matter.
   local idxTar, idxTarP, idxDoP = Lin.IntVec(2), Lin.IntVec(dim), Lin.IntVec(dim)
   idxTar[2], idxTarP[2], idxDoP[2] = yIdxTar, yIdxTar, yIdxTar
   for d=3,dim do idxTarP[d], idxDoP[d] = 1, 1 end   -- Use an arbitrary (valid) index for dirs > 2.

   local xcDo, xcTar, xcP = Lin.Vec(2), Lin.Vec(2), Lin.Vec(dim)

   local yShIndexer, yShPtr = yShift:genIndexer(), yShift:get(1)

   local cellLimTar = {{lo=0., up=0.}, {lo=0., up=0.}}
   local cellLimDo  = {{lo=0., up=0.}, {lo=0., up=0.}}
   local interPts   = {loTar={loDo=nil, upDo=nil}, upTar={loDo=nil, upDo=nil}}
   local dxP        = {}
   for d=1,dim do dxP[d] = 0. end  -- Re-set below.

   local xRange = localRange:selectFirst(1) 
   for xIdx in xRange:rowMajorIter() do
      idxTar[1], idxTarP[1] = xIdx[1], xIdx[1]

      grid:setIndex(idxTarP)
      grid:cellCenter(xcP)
      grid:getDx(dxP)
      for d = 1,2 do xcTar[d], dx[d] = xcP[d], dxP[d] end
      
      cellLimTar = { {lo=grid:cellLowerInDir(1), up=grid:cellUpperInDir(1)},
                     {lo=grid:cellLowerInDir(2), up=grid:cellUpperInDir(2)} }

      yShift:fill(yShIndexer(xIdx), yShPtr)

      local doCellsC = doCells[idxTar[1]][idxTar[2]]   -- Donor cells of current target cell.

      printDebug(string.format("xIdx = %d",xIdx[1]))

      for iC = 1, #doCellsC do

         local idxDo = doCellsC[iC]
         printDebug(string.format("  from idxDo = %d, %d",idxDo[1],idxDo[2]))
         idxDoP[1], idxDoP[2] = idxDo[1], idxDo[2]

         grid:setIndex(idxDoP)
         grid:cellCenter(xcP)
         for d = 1, 2 do xcDo[d] = xcP[d] end

         cellLimDo = { {lo=grid:cellLowerInDir(1), up=grid:cellUpperInDir(1)},
                       {lo=grid:cellLowerInDir(2), up=grid:cellUpperInDir(2)} }

         -- Find the points where y_{j_tar-/+1/2}-yShift intersect the y=y_{j_do-/+1/2} lines.
         interPts = {loTar={loDo=nil, upDo=nil}, upTar={loDo=nil, upDo=nil}}
         local foundAll = true
         for i, _ in pairs(interPts) do
            local tarLim = string.gsub(i, "Tar", "")
            for _, j in ipairs({"loDo","upDo"}) do
               local doLim     = string.gsub(j, "Do", "")
               local yTar, yDo = cellLimTar[2][tarLim], cellLimDo[2][doLim]
               interPts[i][j]  = findIntersect(yTar, yDo, cellLimTar[1], domLim[2])
               if interPts[i][j]==nil then foundAll=false end 
               printDebug(string.format("    O_%s,%s = %g",i,j,interPts[i][j] or -1e19))
            end
         end

         local atUpperCell = idxDo[2]==grid:numCells(2)   -- User in wrapNum

         if foundAll then   -- Scenario sN.
            subCellInt_sN(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
         else

            -- Number of missing intersection points (nil's in
            -- interPts), their indices in interPts, and indices
            -- of found intersection points (non-nil's in interPts).
            local nilNum, nilIdxs, nonNilIdxs = missingInterPt(interPts)
            
            if nilNum == 3 then
               -- si:   y_{j_tar-1/2}-yShift intersects x_{i-1/2}.
               -- sii:  y_{j_tar-1/2}-yShift intersects x_{i+1/2}.
               -- siii: y_{j_tar+1/2}-yShift intersects x_{i-1/2}.
               -- siv:  y_{j_tar+1/2}-yShift intersects x_{i+1/2}.
               if nonNilIdxs[1][1]=="loTar" then    -- Scenario si or sii.
                  subCellInt_siORsii(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
               else   -- Scenario siii or siv.
                  subCellInt_siiiORsiv(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
               end
            elseif nilNum == 1 then
               -- sv:    y_{j_tar+1/2}-yShift doesn't intersect y_{j_do-1/2} & intersects x_{i-1/2}.
               -- svi:   y_{j_tar+1/2}-yShift doesn't intersect y_{j_do-1/2} & intersects x_{i+1/2}.
               -- svii:  y_{j_tar-1/2}-yShift doesn't intersect y_{j_do+1/2} & intersects x_{i-1/2}.
               -- sviii: y_{j_tar-1/2}-yShift doesn't intersect y_{j_do+1/2} & intersects x_{i+1/2}.
               if (nilIdxs[1][1]=="upTar") and (nilIdxs[1][2]=="loDo") then   -- Scenario sv or svi.
                  subCellInt_svORsvi(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
               else   -- Scenario svii or sviii.
                  subCellInt_sviiORsviii(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
               end
            elseif nilNum == 2 then
               -- six:   y_{j_tar-1/2}-yShift crosses y_{j_do-/+1/2} (increasing yShift).
               -- sx:    y_{j_tar-1/2}-yShift crosses y_{j_do-/+1/2} (decreasing yShift).
               -- sxi:   y_{j_tar+1/2}-yShift crosses y_{j_do-/+1/2} (decreasing yShift).
               -- sxii:  y_{j_tar+1/2}-yShift crosses y_{j_do-/+1/2} (increasing yShift).
               -- sxiii: y_{j_tar-1/2}-yShift crosses y_{j_do-1/2} & y_{j_tar+1/2}-yShift crosses y_{j_do+1/2} (increasing yShift).
               -- sxiv:  y_{j_tar-1/2}-yShift crosses y_{j_do-1/2} & y_{j_tar+1/2}-yShift crosses y_{j_do+1/2} (decreasing yShift).
               local nilNumDestLo, nilNumDestUp = 0, 0
               for _, i in ipairs({"loDo","upDo"}) do
                  if interPts.loTar[i] == nil then nilNumDestLo = nilNumDestLo+1 end
                  if interPts.upTar[i] == nil then nilNumDestUp = nilNumDestUp+1 end
               end
               if nilNumDestLo==2 or nilNumDestUp==2 then   -- Scenario six, sx, sxi or sxii.
                  if nilNumDestLo==0 then   -- Scenario six or sx.
                     if interPts.loTar.loDo < interPts.loTar.upDo then    -- Scenario six (similar to scenario si).
                        subCellInt_six(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
                     else   -- Scenario sx.
                        subCellInt_sx(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
                     end
                  else
                     if interPts.upTar.upDo < interPts.upTar.loDo then   -- Scenario sxi.
                        subCellInt_sxi(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
                     else   -- Scenario sxii.
                        subCellInt_sxii(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
                     end
                  end
               else   -- Scenario sxiii or sxiv.
                  subCellInt_sxiiiORsxiv(interPts, xIdx, xcDo, xcTar, cellLimTar, atUpperCell, yShPtr:data(), iC)
               end
            elseif nilNum == 4 then
               -- sxv:  y_{j_tar-1/2}-yShift crosses x_{i-/+1/2}.
               -- sxvi: y_{j_tar+1/2}-yShift crosses x_{i-/+1/2}.
               subCellInt_sxvORsxvi(interPts, xIdx, xcDo, xcTar, cellLimDo, cellLimTar, atUpperCell, yShPtr:data(), iC)
            end
         end
      end

   end
end

return _M
