-- Gkyl ------------------------------------------------------------------------
--
-- Updater to evaluate an analytic function on nodes. 
-- For now it only evaluates on Serendipity nodes.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc            = require "Lib.Alloc"
local SerendipityNodes = require "Lib.SerendipityNodes"
local Lin              = require "Lib.Linalg"
local Proto            = require "Proto"
local Range            = require "Lib.Range"
local UpdaterBase      = require "Updater.Base"

-- System libraries.
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

-- Template for function to map computional space -> physical space.
local compToPhysTempl = xsys.template([[
return function (eta, dx, xc, xOut)
|for i = 1, NDIM do
   xOut[${i}] = 0.5*dx[${i}]*eta[${i}] + xc[${i}]
|end
end
]])

-- Template for function to compute value of function at specified coordinates.
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

-- Projection updater object.
local EvalOnNodes = Proto(UpdaterBase)

function EvalOnNodes:setFunc(func) self._evaluate = func end

function EvalOnNodes:init(tbl)
   EvalOnNodes.super.init(self, tbl) -- Setup base object.

   self._isFirst = true -- For use in advance().

   self._onGrid  = assert(tbl.onGrid, "Updater.EvalOnNodes: Must provide grid object using 'onGrid'")
   self._basis   = assert(tbl.basis, "Updater.EvalOnNodes: Must specify basis functions to use using 'basis'")
   local ev      = assert(tbl.evaluate, "Updater.EvalOnNodes: Must specify function to project using 'evaluate'")
   self:setFunc(ev)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   self._onGhosts = xsys.pickBool(tbl.onGhosts, false)

   local ndim      = self._basis:ndim()
   local polyOrder = self._basis:polyOrder()

   self.nodes    = SerendipityNodes["nodes"..ndim.."xp"..polyOrder]
   self.numNodes = #self.nodes

   -- Construct various functions from template representations.
   self._compToPhys = loadstring(compToPhysTempl {NDIM = ndim} )()

   self.dx  = Lin.Vec(ndim) -- Cell shape.
   self.xc  = Lin.Vec(ndim) -- Cell center.
   self.xmu = Lin.Vec(ndim) -- Coordinate at ordinate.
end

-- Advance method.
function EvalOnNodes:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local qOut = assert(outFld[1], "EvalOnNodes.advance: Must specify an output field")

   local ndim      = grid:ndim()
   local polyOrder = self._basis:polyOrder()
   local numBasis  = self._basis:numBasis()
   local numVal    = qOut:numComponents()/numBasis

   if self._isFirst then
      -- Construct function to evaluate function at specified coordinate.
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
      self.numVal    = self.numVal and self.numVal or numVal
      self.fv        = Lin.Mat(self.numNodes, numVal) -- Function values at ordinates.
   end
   assert(numVal == self.numVal, "EvalOnNodes: created for a scalar/vector field, not a vector/scalar field.")

   local localRangeOut = self._onGhosts and qOut:localExtRange() or qOut:localRange()

   local indexer = qOut:genIndexer()
   local fItr    = qOut:get(1)

   -- Loop, computing projections in each cell.
   for idx in localRangeOut:colMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dx)
      grid:cellCenter(self.xc)

      -- Precompute value of function at each node.
      for i = 1, self.numNodes do
	 self._compToPhys(self.nodes[i], self.dx, self.xc, self.xmu)      -- Compute physical coordinate xmu.
	 self._evalFunc(tCurr, self.xmu, self._evaluate, self.fv[i]) -- Compute function value fv at xmu.
      end

      qOut:fill(indexer(idx), fItr)
      ffiC.nodToMod(self.fv:data(), self.numNodes, numVal, ndim, polyOrder, fItr:data())
   end

   -- Set id of output to id of projection basis.
   qOut:setBasisId(self._basis:id())

   self._isFirst = false
end

return EvalOnNodes
