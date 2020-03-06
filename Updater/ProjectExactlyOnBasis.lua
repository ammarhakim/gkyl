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
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"

-- System libraries.
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

-- scimath: load precomputed ordinates and abscissas for
-- double-exponential integration
local data = require "sci.quad._dblexp_precomputed"
local gmath = require "sci.math".generic

local abscissas = data.abscissas
local weights = data.weigths -- spelling error in scimath!


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
  void projectF(double* f, double* weights, double* basisAtOrdinates, double* fv, int numVal, int numBasis, int numOrd);
]]

-- Projection  updater object.
local ProjectExactlyOnBasis = Proto(UpdaterBase)

function ProjectExactlyOnBasis:setFunc(func)
   self._evaluate = func
end

function ProjectExactlyOnBasis:init(tbl)
   ProjectExactlyOnBasis.super.init(self, tbl) -- Setup base object.

   self._isFirst = true -- For use in advance().

   self._onGrid = assert(tbl.onGrid, "Updater.ProjectExactlyOnBasis: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.ProjectExactlyOnBasis: Must specify basis functions to use using 'basis'")
   local ev = assert(tbl.evaluate, "Updater.ProjectExactlyOnBasis: Must specify function to project using 'evaluate'")
   self:setFunc(ev)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   local ndim = self._basis:ndim()
      
   -- Construct various functions from template representations.
   self._compToPhys = loadstring(compToPhysTempl {NDIM = ndim} )()
end

-- Advance method.
function ProjectExactlyOnBasis:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local qOut = assert(outFld[1], "ProjectExactlyOnBasis.advance: Must specify an output field")

   local ndim = grid:ndim()
   local numOrd  = #abscissas[#abscissas]
   local numBasis = self._basis:numBasis()
   local numVal = qOut:numComponents()/numBasis

   if self._isFirst then
      -- Construct function to evaluate function at specified coorindate.
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
   end

   -- sanity check: ensure number of variables, components and basis functions are consistent
   assert(qOut:numComponents() % numBasis == 0, "ProjectExactlyOnBasis:advance: Incompatible input field")

   local dx = Lin.Vec(ndim) -- cell shape
   local xc = Lin.Vec(ndim) -- cell center
   local fv = Lin.Mat(numOrd, numVal) -- function values at ordinates
   local xmu = Lin.Vec(ndim) -- coordinate at ordinate

   local tId = grid:subGridSharedId() -- local thread ID
   -- object to iterate over only region owned by local SHM thread
   local localRangeDecomp
   if self._projectOnGhosts then
      localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = qOut:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = qOut:localRange(), numSplit = grid:numSharedProcs() }
   end

   local indexer = qOut:genIndexer()
   local fItr = qOut:get(1)

   -- loop, computing projections in each cell
   for idx in localRangeDecomp:colMajorIter(tId) do
      grid:setIndex(idx)
      grid:getDx(dx)
      grid:cellCenter(xc)
   end

   -- set id of output to id of projection basis
   qOut:setBasisId(self._basis:id())

   self._isFirst = false
end

return ProjectExactlyOnBasis
