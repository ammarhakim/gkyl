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
local LinearDecomp     = require "Lib.LinearDecomp"
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

function EvalOnNodes:init(tbl)
   EvalOnNodes.super.init(self, tbl) -- Setup base object.

   self._isFirst  = true -- For use in advance().

   self._onGrid   = assert(tbl.onGrid, "Updater.EvalOnNodes: Must provide grid object using 'onGrid'")
   self._basis    = assert(tbl.basis, "Updater.EvalOnNodes: Must specify basis functions to use using 'basis'")
   self._evaluate = assert(tbl.evaluate, "Updater.EvalOnNodes: Must specify function to project using 'evaluate'")

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   self._projectOnGhosts = xsys.pickBool(tbl.projectOnGhosts, false)

   local ndim      = self._basis:ndim()
   local polyOrder = self._basis:polyOrder()

   self.nodes    = SerendipityNodes["nodes"..ndim.."xp"..polyOrder]
   self.numNodes = #self.nodes

   -- Construct various functions from template representations.
   self._compToPhys = loadstring(compToPhysTempl {NDIM = ndim} )()
end

-- Advance method.
function EvalOnNodes:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local qOut = assert(outFld[1], "EvalOnNodes.advance: Must specify an output field")

   local ndim      = grid:ndim()
   local numBasis  = self._basis:numBasis()
   local polyOrder = self._basis:polyOrder()
   local numVal    = qOut:numComponents()/numBasis

   if self._isFirst then
      -- Construct function to evaluate function at specified coordinate.
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
   end

   -- Sanity check: ensure number of variables, components and basis functions are consistent.
   assert(qOut:numComponents() % numBasis == 0, "EvalOnNodes:advance: Incompatible input field")

   local dx = Lin.Vec(ndim)                  -- Cell shape.
   local xc = Lin.Vec(ndim)                  -- Cell center.
   local fv = Lin.Mat(self.numNodes, numVal) -- Function values at ordinates.
   local xi = Lin.Vec(ndim)                  -- Coordinates at node.

   local tId = grid:subGridSharedId() -- Local thread ID.
   -- Object to iterate over only region owned by local SHM thread.
   local localRangeDecomp
   if self._projectOnGhosts then
      localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = qOut:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = qOut:localRange(), numSplit = grid:numSharedProcs() }
   end

   local indexer = qOut:genIndexer()
   local fItr    = qOut:get(1)

   -- Loop, computing projections in each cell.
   for idx in localRangeDecomp:colMajorIter(tId) do
      grid:setIndex(idx)
      for d = 1, ndim do dx[d] = grid:dx(d) end
      grid:cellCenter(xc)

      -- Precompute value of function at each node.
      for i = 1, self.numNodes do
	 self._compToPhys(self.nodes[i], dx, xc, xi)      -- Compute physical coordinate xi.
	 self._evalFunc(tCurr, xi, self._evaluate, fv[i]) -- Compute function value fv at xi.
      end

      qOut:fill(indexer(idx), fItr)
      ffiC.nodToMod(fv:data(), self.numNodes, numVal, ndim, polyOrder, fItr:data())
   end

   -- Set id of output to id of projection basis.
   qOut:setBasisId(self._basis:id())

   self._isFirst = false
end

return EvalOnNodes
