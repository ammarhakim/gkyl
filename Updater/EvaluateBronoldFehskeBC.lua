local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"

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

   assert(self.grid:ndim() == self.basis:ndim(), "Dimensions of basis and grid must match")
   
   -- Number of quadrature points in each direction
   self.N = 8
   
    -- 1D weights and ordinates
   self.iOrdinates, self.iWeights = GaussQuadRules.ordinates[self.N], GaussQuadRules.weights[self.N]

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
   self.idxP[1] = 1
end

function EvaluateBronoldFehskeBC:_advance(tCurr, inFld, outFld)
   local reflectionFunction = assert(outFld[1], "EvaluateBronoldFehskeBC.advance: Must specify an output field")
   self.vdim = self.ndim - 1
   local Nx = Lin.Vec(self.vdim)
   local dx = Lin.Vec(self.vdim)
   local vc = {}
   local l, u = {}, {}
   for d = 1, self.vdim do
      dx[d] = self.grid:dx(d + 1)
      Nx[d] = self.grid:numCells(d + 1)
      l[d], u[d] = 1, Nx[d]
      vc[d] = {}
      for n = 1, Nx[d] do
	 vc[d][n] = self.grid:lower(d + 1) + (0.5 + (n - 1))*dx[d]
      end
   end
   local Nv = Range.Range(l, u) -- For looping over quadrature nodes.
   local refIndexer = reflectionFunction:genIndexer()
   local fItr = reflectionFunction:get(1)
   for idx in Nv:rowMajorIter() do
      for d = 1, self.vdim do self.idxP[d + 1] = idx[d] end
      reflectionFunction:fill(refIndexer(self.idxP), fItr)
      local vcx = Lin.Vec(self.vdim)
      for d = 1, self.vdim do
	 vcx[d] = vc[d][idx[d]]
      end
      local temp = {}
      for m = 1, self.numBasis do
	 for n = 1, self.numBasis do
	    Rquad = {}
	    for i = 1, self.numOrdinates do
	       local v = Lin.Vec(self.vdim)
	       for d = 1, self.vdim do
	          v[d] = self.ordinates[i][d + 1]
	       end
	       Rquad[i] = self:getR(self:getE(v, vcx, dx), self:getXi(v, vcx, dx))*self.negBasisAtOrdinates[i][n]*self.basisAtOrdinates[i][m]
	    end
	    fItr[self.numBasis*(m - 1) + n] = self:quad(Rquad, -1, 1, self.weights, self.numOrdinates)
	 end
      end
   end
end
      
function EvaluateBronoldFehskeBC:getEta(E, xi) return math.sqrt(1 - ((E - self.electronAffinity)/(self.effectiveMass*E))*(1 - xi^2)) end

function EvaluateBronoldFehskeBC:getE(v, vc, dv)
   local ans = 0
   for d = 1, self.vdim do
      ans = ans + 0.5*self.me*(vc[d] + v[d]*(dv[d]/2))^2/self.elemCharge
   end
   return ans
end

function EvaluateBronoldFehskeBC:getXi(v, vc, dv)
   local ans = 0
   for d = 1, self.vdim do
      ans = ans + (vc[d] + v[d]*(dv[d]/2))^2
   end
   return math.abs(vc[1] + v[1]*(dv[1]/2))/math.sqrt(ans)
end

function EvaluateBronoldFehskeBC:quad(f, a, b, weights, numOrdinates)
   local ans = 0
   for i = 1, numOrdinates do
      ans = ans + (b - a)*(weights[i]*f[i]/2)
   end
   return ans
end

function EvaluateBronoldFehskeBC:getF(E, xi)
   return E^2 - xi^2
end

function EvaluateBronoldFehskeBC:getXic(E)
   if E < (self.electronAffinity/(1 - self.effectiveMass)) then
      return 0.0
   else
      return math.sqrt(1 - ((self.effectiveMass*E)/(E - self.electronAffinity)))
   end
end

function EvaluateBronoldFehskeBC:getT(E, xi)
   local eta = self:getEta(E, xi)
   return (4*self.effectiveMass*math.sqrt(E - self.electronAffinity)*xi*math.sqrt(self.effectiveMass*E)*eta)/(self.effectiveMass*math.sqrt(E - self.electronAffinity)*xi + math.sqrt(self.effectiveMass*E)*eta)^2
end

function EvaluateBronoldFehskeBC:getR(E, xi)
   local xic = self:getXic(E)
   local Tr = self:getT(E, xi)
   local C = self.roughness
   local Tquad = {}
   for i = 1, self.N do
      Tquad[i] = self:getT(E, (1 - xic)*self.iOrdinates[i]/2 + (1 + xic)/2)
   end
   local Tint = self:quad(Tquad, xic, 1.0, self.iWeights, self.N)
   if E < self.electronAffinity then
      return 1.0
   elseif xi > xic then
      return 1 - (Tr/(1 + (C/xi))) - ((C/xi)/(1 + (C/xi)))*Tint
   else
      return 1 - ((C/xi)/(1 + (C/xi)))*Tint
   end
end

return EvaluateBronoldFehskeBC
