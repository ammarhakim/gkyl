-- Gkyl ------------------------------------------------------------------------
--
-- Constant diffusion equation on a rectangular mesh for multimoment models
-- where the moments are store consecutively in every single cell.
--
-- The approach here is to loop over moments and re-use the ConstDiffusion
-- kernels. Probably not the best for memory management, but it's a simpler
-- first step.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local Lin     = require "Lib.Linalg"
local Proto   = require "Lib.Proto"
local ffi     = require "ffi"
local xsys    = require "xsys"
local EqBase  = require "Eq.EqBase"
local ModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"

-- MultimomentConstDiffusion equation on a rectangular mesh.
local MultimomentConstDiffusion = Proto(EqBase)

-- ctor
function MultimomentConstDiffusion:init(tbl)

   self.nMoments = assert(tbl.numMoments, "Eq.MultimomentConstDiffusion: Must indicate number of moments in 'numMoments'.")

   self._basis = assert(tbl.basis,
      "Eq.MultimomentConstDiffusion: Must specify basis functions to use using 'basis'")

   local nm, dim, pOrder = self._basis:id(), self._basis:ndim(), self._basis:polyOrder()
   self.numB = self._basis:numBasis()

   -- Read the directions in which to apply diffusion. 
   local diffDirsIn = tbl.diffusiveDirs
   if diffDirsIn then
      assert(#diffDirsIn<=dim, "Eq.MultimomentConstDiffusion: 'diffusiveDirs' cannot have more entries than the simulation's dimensions.")
      self.diffDirs = diffDirsIn
   else
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, dim do self.diffDirs[d] = d end
   end
   
   -- Read diffusion coefficient (or vector).
   local nuIn     = assert(tbl.coefficient,
                           "Eq.MultimomentConstDiffusion: must specify diffusion coefficient (or vector) using 'coefficient' ")
   local nuInType = type(nuIn)
   local isVarCoeff = false
   if (nuInType == "number") then
      -- Set the diffusion coefficient to the same amplitude in all directions.
      self._nu = Lin.Vec(dim)
      for d = 1, dim do self._nu[d] = nuIn end
   elseif (nuInType == "table") then
      if #nuIn > 0 then   -- Must be a table of numbers.
         if diffDirsIn then
            assert(#nuIn==#diffDirsIn, "Eq.MultimomentConstDiffusion: 'coefficient' table must have the same number of entries as 'diffusiveDirs'.")
         else
            assert(#nuIn==dim, "Eq.MultimomentConstDiffusion: 'coefficient' table must have the same number of entries as the simulation's dimensions if 'diffusiveDirs' is not given.")
         end
         self._nu = Lin.Vec(dim)
         for d = 1, dim do self._nu[d] = 0.0 end
         for d = 1, #self.diffDirs do self._nu[self.diffDirs[d]] = nuIn[d] end
      else

         -- The coefficient is a CartField with the spatially varying coefficient.
         isVarCoeff = true
         self._nu   = nuIn
      end
   else
      assert(false, "Eq.MultimomentConstDiffusion: 'coefficient' must be a number or a table.")
   end

   local diffOrder
   if tbl.order then
      diffOrder = tbl.order
      assert(not (diffOrder > 4 and pOrder < 2), "Eq.ConstDiffusion: grad^6 requires polyOrder > 1.")
   else
      diffOrder = 2
   end

   local applyPositivity = xsys.pickBool(tbl.positivity, false)   -- Positivity preserving option.

   -- Store pointers to C kernels implementing volume and surface terms.
   self._volTerm           = ModDecl.selectVol(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff)
   self._surfTerms         = ModDecl.selectSurf(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff, applyPositivity)
   self._boundarySurfTerms = ModDecl.selectBoundarySurf(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff, applyPositivity)
   
   if isVarCoeff then

      self._nuPtrl = self._nu:get(1)
      self._nuPtrr = self._nu:get(1)
      self._nuIdxr = self._nu:genIndexer()

      self.volTermFunc = function(w, dx, idx, q, out)
         return MultimomentConstDiffusion["volTermVarInSpace"](self, w, dx, idx, q, out)
      end
      self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return MultimomentConstDiffusion["surfTermVarInSpace"](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
      self.boundarySurfTermFunc = function(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return MultimomentConstDiffusion["boundarySurfTermVarInSpace"](self, dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
   else
      self.volTermFunc = function(w, dx, idx, q, out)
         return MultimomentConstDiffusion["volTermConstInSpace"](self, w, dx, idx, q, out)
      end
      self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return MultimomentConstDiffusion["surfTermConstInSpace"](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
      self.boundarySurfTermFunc = function(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return MultimomentConstDiffusion["boundarySurfTermConstInSpace"](self, dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
   end
end

-- Methods.
function MultimomentConstDiffusion:numEquations() return 1 end
function MultimomentConstDiffusion:numWaves() return 1 end
function MultimomentConstDiffusion:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- Flux in direction dir.
function MultimomentConstDiffusion:flux(dir, qIn, fOut)
   assert(false, "MultimomentConstDiffusion:flux: NYI!")
end

-- Riemann problem for MultimomentConstDiffusion equation.
function MultimomentConstDiffusion:rp(dir, delta, ql, qr, waves, s)
   assert(false, "MultimomentConstDiffusion:rp: NYI!")
end

-- Compute q-fluctuations.
function MultimomentConstDiffusion:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "MultimomentConstDiffusion:qFluctuations: NYI!")
end

-- Maximum wave speed.
function MultimomentConstDiffusion:maxSpeed(dir, w, dx, q)
   assert(false, "MultimomentConstDiffusion:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function MultimomentConstDiffusion:volTerm(w, dx, idx, q, out)
   return self.volTermFunc(w, dx, idx, q, out)
end
function MultimomentConstDiffusion:volTermConstInSpace(w, dx, idx, q, out)
   local cflFreq = 0.
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      cflFreq = math.max(cflFreq,self._volTerm(w:data(), dx:data(), self._nu:data(), q:data(), out:data()))
   end
   return cflFreq
end
function MultimomentConstDiffusion:volTermVarInSpace(w, dx, idx, q, out)
   self._nu:fill(self._nuIdxr(idx), self._nuPtrl)
   local cflFreq = 0.
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      cflFreq = math.max(cflFreq,self._volTerm(w:data(), dx:data(), self._nuPtrl:data(), q:data()+momOff, out:data())+momOff)
   end
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function MultimomentConstDiffusion:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
end
function MultimomentConstDiffusion:surfTermConstInSpace(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nu:data(), ql:data()+momOff, qr:data()+momOff, outl:data()+momOff, outr:data()+momOff)
   end
   return 0
end
function MultimomentConstDiffusion:surfTermVarInSpace(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._nu:fill(self._nuIdxr(idxl), self._nuPtrl)
   self._nu:fill(self._nuIdxr(idxr), self._nuPtrr)
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuPtrl:data(), self._nuPtrr:data(), ql:data()+momOff, qr:data()+momOff, outl:data()+momOff, outr:data()+momOff)
   end
   return 0
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function MultimomentConstDiffusion:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   return self.boundarySurfTermFunc(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
end
function MultimomentConstDiffusion:boundarySurfTermConstInSpace(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nu:data(), ql:data()+momOff, qr:data()+momOff, outl:data()+momOff, outr:data()+momOff)
   end
   return 0
end
function MultimomentConstDiffusion:boundarySurfTermVarInSpace(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._nu:fill(self._nuIdxr(idxl), self._nuPtrl)
   self._nu:fill(self._nuIdxr(idxr), self._nuPtrr)
   for mI = 1, self.nMoments do
      local momOff = (mI-1)*self.numB
      self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), self._nuPtrl:data(), self._nuPtrr:data(), ql:data()+momOff, qr:data()+momOff, outl:data()+momOff, outr:data()+momOff)
   end
   return 0
end

function MultimomentConstDiffusion:setAuxFields(auxFields) end


return MultimomentConstDiffusion
