-- Gkyl ------------------------------------------------------------------------
--
-- Updater to project function on basis functions. Uses Gaussian
-- quadrature. The projection is exact if the function being projected
-- is a polynomial of order less than the basis function polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc          = require "Lib.Alloc"
local Lin            = require "Lib.Linalg"
local Proto          = require "Proto"
local UpdaterBase    = require "Updater.Base"
local ZeroArray      = require "DataStruct.ZeroArray"

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
// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD = 0, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD, // Gauss-Lobatto quadrature
};

// Object type
typedef struct gkyl_proj_on_basis gkyl_proj_on_basis;

// input packaged as a struct
struct gkyl_proj_on_basis_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project
  const struct gkyl_basis *basis; // basis functions

  enum gkyl_quad_type qtype; // quadrature to use

  int num_quad; // number of quadrature points
  int num_ret_vals; // number of return values in eval function
  evalf_t eval; // function to project
  void *ctx; // function context
};

/**
 * Create new updater to project function on basis functions on a
 * grid. Free using gkyl_proj_on_basis_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_quad Number of quadrature nodes
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
struct gkyl_proj_on_basis *gkyl_proj_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, int num_quad, int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Create new updater to project function on basis functions on a
 * grid. Free using gkyl_proj_on_basis_release method.
 *
 * @param inp Input parameters
 * @return New updater pointer.
 */
struct gkyl_proj_on_basis* gkyl_proj_on_basis_inew(const struct gkyl_proj_on_basis_inp *inp);

/**
 * Perform the quadrature in the proj_on_basis procedure.
 * Intended for systems that can't perform the whole procedure in _advance.
 *
 * @param up Project on basis updater.
 * @param fun_at_ords Function evaluated at ordinates in one cell.
 * @param f Output projected function in one cell.
 */
void gkyl_proj_on_basis_quad(const struct gkyl_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f);

/**
 * Return the total number of quadrature points/ordinates.
 *
 * @param up Project on basis updater.
 * @return Number of ordinates.
 */
int gkyl_proj_on_basis_get_tot_quad(const struct gkyl_proj_on_basis *up);

/**
 * Get the coordinates of a given ordinate.
 *
 * @param up Project on basis updater.
 * @param node Index indicate the desired node.
 * @return Node coordinates.
 */
double* gkyl_proj_on_basis_fetch_ordinate(const struct gkyl_proj_on_basis *up, long node);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_on_basis_release(struct gkyl_proj_on_basis* pob);
]]

-- Projection  updater object.
local ProjectOnBasis = Proto(UpdaterBase)

function ProjectOnBasis:setFunc(func) self._evaluate = func end

function ProjectOnBasis:init(tbl)
   ProjectOnBasis.super.init(self, tbl) -- Setup base object.

   self._isFirst = true -- For use in advance().

   self._onGrid = assert(tbl.onGrid, "Updater.ProjectOnBasis: Must provide grid object using 'onGrid'")
   self._basis  = assert(tbl.basis, "Updater.ProjectOnBasis: Must specify basis functions to use using 'basis'")
   local ev     = assert(tbl.evaluate, "Updater.ProjectOnBasis: Must specify function to project using 'evaluate'")
   self:setFunc(ev)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   self.quadType = tbl.rule and tbl.rule or "gauss_legendre" -- Number of quadrature points in each direction
   -- Number of quadrature points in each direction.
   self.numQuad  = tbl.numQuad and tbl.numQuad or (
     self.quadType=="gauss_legendre" and self._basis:polyOrder()+1 or self._basis:polyOrder()+2 
   )

   self._onGhosts = xsys.pickBool(tbl.onGhosts, false)

   local ndim = self._basis:ndim()

   -- Construct various functions from template representations.
   self._compToPhys = loadstring(compToPhysTempl {NDIM = ndim} )()

   self.dx  = Lin.Vec(ndim) -- Cell shape.
   self.xc  = Lin.Vec(ndim) -- Cell center.
   self.xmu = Lin.Vec(ndim) -- Coordinate at ordinate.
end

-- Advance method.
function ProjectOnBasis:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local qOut = assert(outFld[1], "ProjectOnBasis.advance: Must specify an output field")

   local numBasis = self._basis:numBasis()
   local numVal   = qOut:numComponents()/numBasis

   if self._isFirst then
      -- Construct function to evaluate function at specified coordinate.
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
      self.numVal    = self.numVal and self.numVal or numVal
      local projInp  = ffi.new("struct gkyl_proj_on_basis_inp")
      projInp.grid         = self._onGrid._zero
      projInp.basis        = self._basis._zero
      projInp.qtype        = self.quadType=="gauss_legendre" and 0 or 1
      projInp.num_quad     = self.numQuad
      projInp.num_ret_vals = self.numVal
      projInp.eval         = nil
      projInp.ctx          = nil
      self._zero = ffi.gc(ffiC.gkyl_proj_on_basis_inew(projInp),
                          ffiC.gkyl_proj_on_basis_release)
      self.numNodes  = ffiC.gkyl_proj_on_basis_get_tot_quad(self._zero)
      self.fv        = ZeroArray.Array(ZeroArray.double, numVal, self.numNodes, 0)
   end
   assert(numVal == self.numVal, "ProjectOnBasis: created for a scalar/vector field, not a vector/scalar field.")

   -- Sanity check: ensure number of variables, components and basis functions are consistent.
   assert(qOut:numComponents() % numBasis == 0, "ProjectOnBasis:advance: Incompatible input field")

   -- Object to iterate over only region owned by local SHM thread.
   local localRangeOut = self._onGhosts and qOut:localExtRange() or qOut:localRange()

   local indexer = qOut:genIndexer()
   local fItr    = qOut:get(1)

   -- Loop, computing projections in each cell.
   for idx in localRangeOut:rowMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dx)
      grid:cellCenter(self.xc)

      -- Precompute value of function at each ordinate.
      for i = 1, self.numNodes do
         self._compToPhys(ffiC.gkyl_proj_on_basis_fetch_ordinate(self._zero,i-1), self.dx, self.xc, self.xmu) -- Compute coordinate.
         self._evalFunc(tCurr, self.xmu, self._evaluate, self.fv:fetch(i-1))
      end
      
      qOut:fill(indexer(idx), fItr)
      ffiC.gkyl_proj_on_basis_quad(self._zero, self.fv, fItr:data())
   end

   -- Set id of output to id of projection basis.
   qOut:setBasisId(self._basis:id())

   self._isFirst = false
end

return ProjectOnBasis
