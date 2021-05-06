-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the Bronold & Fehske reflection function for
-- each cell in velocity space. Uses Gaussian quadrature to integrate
-- over grid dimensions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

local ffi  = require "ffi"
local ffiC = ffi.C

ffi.cdef[[
  void reflectionLoop(double* fItr, double* ordinates, double* weights, double* iWeights, double* iOrdinates, double* negBasis, double* posBasis, double* vc, double* dv, double C, double elemCharge, double me, double electronAffinity, double effectiveMass, int vdim, int cdim, int numOrd, int N, int numBasis, int ignore, double* Rquad, double* v, double* Tquad);
]]

-- Updater object.
local EvaluateBronoldFehskeBC = Proto(UpdaterBase)

function EvaluateBronoldFehskeBC:init(tbl)
   EvaluateBronoldFehskeBC.super.init(self, tbl)
   
   self.grid = assert(tbl.onGrid, "Updater.EvaluateBronoldFehskeBC: Must provide grid object using 'onGrid'")
   self.basis  = assert(tbl.basis, "Updater.EvaluateBronoldFehskeBC: Must specify basis functions to use using 'basis'")
   self.electronAffinity = assert(tbl.electronAffinity, "Updater.EvaluateBronoldFehskeBC: Must specify electronAffinity")
   self.effectiveMass = assert(tbl.effectiveMass, "Updater.EvaluateBronoldFehskeBC: Must specify effectiveMass")
   self.elemCharge = assert(tbl.elemCharge, "Updater.EvaluateBronoldFehskeBC: Must specify elemCharge")
   self.me = assert(tbl.me, "Updater.EvaluateBronoldFehskeBC: Must specify electronMass")
   self.roughness = 2.0
   self.vdim = tbl.vdim
   self.cdim = tbl.cdim
   
   -- Option to ignore the negative or positive part of velocity space
   if tbl.ignore == "neg" then self.ignore = -1 elseif tbl.ignore == "pos" then self.ignore = 1 else self.ignore = 0 end

   assert(self.grid:ndim() == self.basis:ndim(), "Dimensions of basis and grid must match")
   
   -- Number of quadrature points in each direction
   self.N = 8
   
   -- 1D weights and ordinates
   self.iOrdinates, self.iWeights = GaussQuadRules.ordinates[self.N], GaussQuadRules.weights[self.N]
   self.iordinates, self.iweights = Lin.Vec(self.N), Lin.Vec(self.N)
   self.ndim = self.basis:ndim()
   local l, u = {}, {}
   for d = 1, self.ndim do l[d], u[d] = 1, self.N end
   local quadRange = Range.Range(l, u) -- For looping over quadrature nodes.

   self.numOrdinates = quadRange:volume() -- Number of ordinates.
   
   -- Construct weights and ordinates for integration in multiple dimensions.
   self.ordinates = Lin.Mat(self.numOrdinates, self.ndim)
   local negOrdinates = Lin.Mat(self.numOrdinates, self.ndim)
   self.weights = Lin.Vec(self.numOrdinates)
   local nodeNum = 1
   for idx in quadRange:rowMajorIter() do
      self.weights[nodeNum] = 1.0
      for d = 1, self.ndim do
	 self.weights[nodeNum] = self.weights[nodeNum]*self.iWeights[idx[d]]
	 self.iweights[idx[d]] = self.iWeights[idx[d]]
	 self.iordinates[idx[d]] = self.iOrdinates[idx[d]]
	 self.ordinates[nodeNum][d] = self.iOrdinates[idx[d]]
	 if d < 3 then
	    negOrdinates[nodeNum][d] = -self.iOrdinates[idx[d]]
	 else
	    negOrdinates[nodeNum][d] = self.iOrdinates[idx[d]]
	 end
      end
      nodeNum = nodeNum + 1
   end

   self.numBasis = self.basis:numBasis()
   self.basisAtOrdinates = Lin.Mat(self.numOrdinates, self.numBasis)
   self.negBasisAtOrdinates = Lin.Mat(self.numOrdinates, self.numBasis)
   
   -- Pre-compute values of basis functions at quadrature nodes.
   if self.numBasis > 1 then
      for n = 1, self.numOrdinates do
	 self.basis:evalBasis(self.ordinates[n], self.basisAtOrdinates[n])
	 self.basis:evalBasis(negOrdinates[n], self.negBasisAtOrdinates[n])
      end
   else
      for n = 1, self.numOrdinates do
	 self.basisAtOrdinates[n][1] = 1.0/2^self.grid:ndim()
      end
   end
   self.idxP = Lin.IntVec(self.ndim)
   for d = 1, self.cdim do
      self.idxP[d] = 1
   end
end

function EvaluateBronoldFehskeBC:_advance(tCurr, inFld, outFld)
   local reflectionFunction = assert(outFld[1], "EvaluateBronoldFehskeBC.advance: Must specify an output field")
   local Nx = Lin.Vec(self.vdim)
   local l, u = {}, {}
   for d = 1, self.vdim do
      Nx[d] = self.grid:numCells(d + self.cdim)
      l[d], u[d] = 1, Nx[d]
   end
   local Nv = Range.Range(l, u) -- For looping over quadrature nodes.
   local refIndexer = reflectionFunction:genIndexer()
   local fItr = reflectionFunction:get(1)

   local xc = Lin.Vec(self.ndim)
   local dx = Lin.Vec(self.ndim)
   local vcx = Lin.Vec(self.vdim)
   local dv = Lin.Vec(self.vdim)
   local Rquad = Lin.Vec(self.numOrdinates)
   local v = Lin.Vec(self.vdim)
   local Tquad = Lin.Vec(self.N)
   
   -- Loop over velocity space
   for idx in Nv:rowMajorIter() do
      for d = 1, self.vdim do self.idxP[d + self.cdim] = idx[d] end
      self.grid:setIndex(self.idxP)
      self.grid:cellCenter(xc)
      self.grid:getDx(dx)
      for d = 1, self.vdim do
         vcx[d] = xc[d + self.cdim]
         dv[d] = dx[d + self.cdim]
      end
      reflectionFunction:fill(refIndexer(self.idxP), fItr)
      ffiC.reflectionLoop(fItr:data(), self.ordinates:data(), self.weights:data(), self.iweights:data(), self.iordinates:data(), self.negBasisAtOrdinates:data(), self.basisAtOrdinates:data(), vcx:data(), dv:data(), self.roughness, self.elemCharge, self.me, self.electronAffinity, self.effectiveMass, self.vdim, self.cdim, self.numOrdinates, self.N, self.numBasis, self.ignore, Rquad:data(), v:data(), Tquad:data())
   end
end

return EvaluateBronoldFehskeBC
