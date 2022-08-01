-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create a bi-Maxwellian distribution from the number density,
-- flow speed and parallel and perpendicular temperatures 
-- and project it on basis functions. Uses Gaussian
-- quadrature.
--
-- Note:
--  - Primarily intended for gyrokinetics for now.
--  - Unlike MaxwellianOnBasis.lua this does not have a "Lua implementation",
--    i.e. the loops are unrolled and the FLOPs are all in C++.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local LinearDecomp = require "Lib.LinearDecomp"
local Lin          = require "Lib.Linalg"
local Proto        = require "Proto"
local Range        = require "Lib.Range"
local Time         = require "Lib.Time"
local UpdaterBase  = require "Updater.Base"
local ModDecl      = require "Updater.maxwellianOnBasisData.BiMaxwellianOnBasisModDecl" 

-- System libraries.
local xsys = require "xsys"

-- Inherit the base Updater from UpdaterBase updater object.
local BiMaxwellianOnBasis = Proto(UpdaterBase)

function BiMaxwellianOnBasis:init(tbl)
   BiMaxwellianOnBasis.super.init(self, tbl) -- setup base object

   self.phaseGrid  = assert(tbl.onGrid,
                            "Updater.BiMaxwellianOnBasis: Must provide phase space grid object 'onGrid'")
   self.confGrid   = assert(tbl.confGrid,
                            "Updater.BiMaxwellianOnBasis: Must provide configuration space grid object 'confGrid'")
   self.confBasis  = assert(tbl.confBasis,
                            "Updater.BiMaxwellianOnBasis: Must provide configuration space basis object 'confBasis'")
   self.phaseBasis = assert(tbl.phaseBasis,
                            "Updater.BiMaxwellianOnBasis: Must provide phase space basis object 'phaseBasis'")

   -- Number of quadrature points in each direction
   local numQuad1D = self.confBasis:polyOrder() + 1

   self.quadType = "Gauss"

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   -- Mass and magnetic field amplitude are needed to project a gyrokinetic Maxwellian.
   self.mass = assert(tbl.mass,
                      "Updater.BiMaxwellianOnBasis: Must provide the species' mass in 'mass'.")

   self._pDim, self._cDim = self.phaseBasis:ndim(), self.confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   self.numConfOrds = numQuad1D^self._cDim 

   self.uFlowOrd = Lin.Vec(self.numConfOrds)
   self.vtparSqOrd  = Lin.Vec(self.numConfOrds)
   self.vtperpSqOrd = Lin.Vec(self.numConfOrds)
   self.normFacOrd  = Lin.Vec(self.numConfOrds)
   self.bmagOrd     = Lin.Vec(self.numConfOrds)

   local quadFuncs = ModDecl.selectQuad(self.confBasis:id(), self._cDim, self._vDim, self.confBasis:polyOrder(), self.quadType, self.isGK)
   self._evAtConfOrds = quadFuncs[1]
   self._phaseQuad    = quadFuncs[2]

   -- Cell index, center, and dx.
   self.idx = Lin.IntVec(self._pDim)
   self.xc  = Lin.Vec(self._pDim)
   self.dx  = Lin.Vec(self._pDim)
end

function BiMaxwellianOnBasis:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs.
   local nIn        = assert(inFld[1], "BiMaxwellianOnBasis: Must specify density 'inFld[1]'")
   local uFlowIn    = assert(inFld[2], "BiMaxwellianOnBasis: Must specify drift speed 'inFld[2]'")
   local vtparSqIn  = assert(inFld[3], "BiMaxwellianOnBasis: Must specify parallel thermal velocity squared 'inFld[3]'")
   local vtperpSqIn = assert(inFld[4], "BiMaxwellianOnBasis: Must specify perpendicular thermal velocity squared 'inFld[4]'")
   local fOut       = assert(outFld[1], "BiMaxwellianOnBasis: Must specify an output field 'outFld[1]'")

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local nItr, uFlowItr, vtparSqItr, vtperpSqItr = nIn:get(1), uFlowIn:get(1), vtparSqIn:get(1), vtperpSqIn:get(1)
   local fItr = fOut:get(1)

   local bmagIn, bmagItr
   bmagIn  = assert(inFld[5], "BiMaxwellianOnBasis: Must specify magnetic field amplitude 'inFld[5]'")
   bmagItr = bmagIn:get(1)

   local confIndexer  = nIn:genIndexer()
   local phaseIndexer = fOut:genIndexer()
   local phaseRange
   if self.onGhosts then
      phaseRange = fOut:localExtRange()
   else
      phaseRange = fOut:localRange()
   end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = self.phaseGrid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = self.phaseGrid:subGridSharedId()   -- Local thread ID.

   -- The configuration space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      nIn:fill(confIndexer(cIdx), nItr)
      uFlowIn:fill(confIndexer(cIdx), uFlowItr)
      vtparSqIn:fill(confIndexer(cIdx), vtparSqItr)
      vtperpSqIn:fill(confIndexer(cIdx), vtperpSqItr)
      bmagIn:fill(confIndexer(cIdx), bmagItr)

      self._evAtConfOrds(nItr:data(), uFlowItr:data(), vtparSqItr:data(), vtperpSqItr:data(), bmagItr:data(),
                         self.uFlowOrd:data(), self.vtparSqOrd:data(), self.vtperpSqOrd:data(), self.normFacOrd:data(), self.bmagOrd:data())

      -- The velocity space loop
      for vIdx in velRange:rowMajorIter() do
         cIdx:copyInto(self.idx)
         for d = 1, vDim do self.idx[cDim+d] = vIdx[d] end

         fOut:fill(phaseIndexer(self.idx), fItr)
   
         -- Get cell shape, cell center coordinates
         self.phaseGrid:setIndex(self.idx)
         self.phaseGrid:getDx(self.dx)
         self.phaseGrid:cellCenter(self.xc)
   
         self._phaseQuad(self.uFlowOrd:data(), self.vtparSqOrd:data(), self.vtperpSqOrd:data(), self.normFacOrd:data(),
                         self.bmagOrd:data(), self.mass,
                         self.xc:data(), self.dx:data(), fItr:data())

      end
   end
   -- Set id of output to id of projection basis.
   fOut:setBasisId(self.phaseBasis:id())
end

return BiMaxwellianOnBasis
