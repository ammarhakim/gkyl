-- Gkyl ------------------------------------------------------------------------
--
-- Project the square root (to a power) of a CartField (or its reciprocal) onto
-- the basis using quadrature. So can compute (sqrt(f))^q where q can be negative.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local GaussQuadRules     = require "Lib.GaussQuadRules"
local LinearDecomp       = require "Lib.LinearDecomp"
local Lin                = require "Lib.Linalg"
local Proto              = require "Proto"
local Range              = require "Lib.Range"
local Time               = require "Lib.Time"
local UpdaterBase        = require "Updater.Base"
local SqrtOnBasisModDecl = require "Updater.sqrtOnBasisData.sqrt_on_basis_mod_decl"

-- System libraries.
local xsys = require "xsys"

-- Inherit the base Updater from UpdaterBase updater object.
local SqrtOnBasis = Proto(UpdaterBase)

function SqrtOnBasis:init(tbl)
   SqrtOnBasis.super.init(self, tbl) -- setup base object

   self.grid  = assert(tbl.onGrid, "Updater.SqrtOnBasis: Must provide configuration space grid object 'grid'.")
   self.basis = assert(tbl.basis, "Updater.SqrtOnBasis: Must provide configuration space basis object 'basis'.")
   local calcReciprocal = tbl.reciprocal or false

   -- Number of quadrature points in each direction
   local numQuad1D = self.basis:polyOrder() + 1

   self.quadType = "gauss"

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self.dim = self.basis:ndim()

   self._sqrtOBKer = SqrtOnBasisModDecl.selectQuad(self.basis:id(), self.dim, self.basis:polyOrder(), self.quadType)
end

function SqrtOnBasis:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs.
   local fIn  = assert(inFld[1], "SqrtOnBasis: Must specify an input field 'inFld[1]'.")
   local fOut = assert(outFld[1], "SqrtOnBasis: Must specify an output field 'outFld[1]'.")

   local exponent = inFld[2] or 1.0

   local dim, grid = self.dim, self.grid

   local fInItr, fOutItr = fIn:get(1), fOut:get(1)

   local indexer = fIn:genIndexer()
   local fOutRangeDecomp
   if self.onGhosts then
      fOutRangeDecomp = LinearDecomp.LinearDecompRange {
         range = fOut:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      fOutRangeDecomp = LinearDecomp.LinearDecompRange {
         range = fOut:localRange(), numSplit = grid:numSharedProcs() }
   end

   local tId = grid:subGridSharedId()   -- Local thread ID.
   for idx in fOutRangeDecomp:rowMajorIter(tId) do
      fIn:fill(indexer(idx), fInItr)
      fOut:fill(indexer(idx), fOutItr)

      self._sqrtOBKer(exponent, fInItr:data(),fOutItr:data())
   end
end

return SqrtOnBasis
