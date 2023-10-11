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
local Lin              = require "Lib.Linalg"
local Proto            = require "Proto"
local UpdaterBase      = require "Updater.Base"
local ZeroArray        = require "DataStruct.ZeroArray"

-- System libraries.
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
require "Lib.ZeroUtil"

-- Template for function to map computional space -> physical space.
local compToPhysTempl = xsys.template([[
return function (eta, dx, xc, xOut)
|for i = 1, NDIM do
   xOut[${i}] = 0.5*dx[${i}]*eta[${i-1}] + xc[${i}]
|end
end
]])

-- Template for function to compute value of function at specified coordinates.
local evalFuncTempl = xsys.template([[
return function (tCurr, xn, func, fout)
|for i = 1, M-1 do
   fout[${i-1}], 
|end
   fout[${M-1}] = func(tCurr, xn)
end
]])

ffi.cdef[[
// Object type
typedef struct gkyl_eval_on_nodes gkyl_eval_on_nodes;

/**
 * Create new updater to compute function on nodes and calculate its
 * expansion on basis functions. Free using gkyl_eval_on_nodes_release
 * method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_eval_on_nodes* gkyl_eval_on_nodes_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute evaluation on nodes and corresponding expansion
 * coefficients. The update_rng MUST be a sub-range of the range on
 * which the array is defined. That is, it must be either the same
 * range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Eval on nodes updater to run
 * @param tm Time at which eval must be computed
 * @param update_rng Range on which to run eval.
 * @param out Output array
 */
void gkyl_eval_on_nodes_advance(const gkyl_eval_on_nodes *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Perform the nodal to modal transformation.
 *
 * @param up Project on basis updater.
 * @param fun_at_nodes Function evaluated at nodes in one cell.
 * @param f Modal coefficients of the function in one cell.
 */
void gkyl_eval_on_nodes_nod2mod(const gkyl_eval_on_nodes *up, const struct gkyl_array *fun_at_nodes, double *f);

/**
 * Get the coordinates of a given node.
 *
 * @param up Project on basis updater.
 * @param node Index indicate the desired node.
 * @return Node coordinates.
 */
double* gkyl_eval_on_nodes_fetch_node(const gkyl_eval_on_nodes *up, long node);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_eval_on_nodes_release(gkyl_eval_on_nodes *up);
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
   self.globalUpdateRange = tbl.globalUpdateRange


   local ndim = self._basis:ndim()

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

   local numBasis = self._basis:numBasis()
   local numVal   = qOut:numComponents()/numBasis

   if self._isFirst then
      -- Construct function to evaluate function at specified coordinate.
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
      self.numVal    = self.numVal and self.numVal or numVal
      self._zero     = ffi.gc(ffiC.gkyl_eval_on_nodes_new(self._onGrid._zero, self._basis._zero, self.numVal, nil, nil),
                              ffiC.gkyl_eval_on_nodes_release)
      self.numNodes  = self._basis:numBasis()
      self.fv        = ZeroArray.Array(ZeroArray.double, numVal, self.numNodes, 0)
   end
   assert(numVal == self.numVal, "EvalOnNodes: created for a scalar/vector field, not a vector/scalar field.")


   local localRangeOut
   if self.globalUpdateRange then
      localRangeOut = qOut:localExtRange():intersect(self.globalUpdateRange)
   else
      localRangeOut = self._onGhosts and qOut:localExtRange() or qOut:localRange()
   end


   local indexer = qOut:genIndexer()
   local fItr    = qOut:get(1)

   -- Loop, computing projections in each cell.
   for idx in localRangeOut:colMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dx)
      grid:cellCenter(self.xc)

      -- Precompute value of function at each node.
      for i = 1, self.numNodes do
	 self._compToPhys(ffiC.gkyl_eval_on_nodes_fetch_node(self._zero,i-1), self.dx, self.xc, self.xmu)      -- Compute physical coordinate xmu.
	 self._evalFunc(tCurr, self.xmu, self._evaluate, self.fv:fetch(i-1)) -- Compute function value fv at xmu.
      end

      qOut:fill(indexer(idx), fItr)
      ffiC.gkyl_eval_on_nodes_nod2mod(self._zero, self.fv, fItr:data())
   end

   -- Set id of output to id of projection basis.
   qOut:setBasisId(self._basis:id())

   self._isFirst = false
end

return EvalOnNodes
