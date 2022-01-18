-- Gkyl ------------------------------------------------------------------------
--
-- Constant diffusion equation on a rectangular mesh.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local Lin    = require "Lib.Linalg"
local Proto  = require "Lib.Proto"
local ffi    = require "ffi"
local xsys   = require "xsys"
local EqBase = require "Eq.EqBase"
local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"

-- ConstDiffusion equation on a rectangular mesh.
local ConstDiffusion = Proto(EqBase)

-- ctor
function ConstDiffusion:init(tbl)

   self._basis = assert(tbl.basis,
      "Eq.constDiffusion: Must specify basis functions to use using 'basis'")

   local nm, dim, pOrder = self._basis:id(), self._basis:ndim(), self._basis:polyOrder()

   -- Read the directions in which to apply diffusion. 
   local diffDirsIn = tbl.diffusiveDirs
   if diffDirsIn then
      assert(#diffDirsIn<=dim, "Eq.constDiffusion: 'diffusiveDirs' cannot have more entries than the simulation's dimensions.")
      self.diffDirs = diffDirsIn
   else
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, dim do self.diffDirs[d] = d end
   end
   
   -- Read diffusion coefficient (or vector).
   local nuIn     = assert(tbl.coefficient,
                           "Eq.constDiffusion: must specify diffusion coefficient (or vector) using 'coefficient' ")
   local nuInType = type(nuIn)
   local isVarCoeff = false
   if (nuInType == "number") then
      -- Set the diffusion coefficient to the same amplitude in all directions.
      self._nu = Lin.Vec(dim)
      for d = 1, dim do self._nu[d] = nuIn end
   elseif (nuInType == "table") then
      if #nuIn > 0 then   -- Must be a table of numbers.
         if diffDirsIn then
            assert(#nuIn==#diffDirsIn, "Eq.constDiffusion: 'coefficient' table must have the same number of entries as 'diffusiveDirs'.")
         else
            assert(#nuIn==dim, "Eq.constDiffusion: 'coefficient' table must have the same number of entries as the simulation's dimensions if 'diffusiveDirs' is not given.")
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
      assert(false, "Eq.constDiffusion: 'coefficient' must be a number or a table.")
   end

   local diffOrder
   if tbl.order then
      diffOrder = tbl.order
      assert(not (diffOrder > 4 and pOrder < 2), "Eq.constDiffusion: grad^6 requires polyOrder > 1.")
   else
      diffOrder = 2
   end

   local applyPositivity = xsys.pickBool(tbl.positivity, false)   -- Positivity preserving option.

   -- Store pointers to C kernels implementing volume and surface terms.
   self._volTerm           = ConstDiffusionModDecl.selectVol(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff)
   self._surfTerms         = ConstDiffusionModDecl.selectSurf(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff, applyPositivity)
   self._boundarySurfTerms = ConstDiffusionModDecl.selectBoundarySurf(nm, dim, pOrder, self.diffDirs, diffOrder, isVarCoeff, applyPositivity)
   
   if isVarCoeff then

      self._nuPtrl = self._nu:get(1)
      self._nuPtrr = self._nu:get(1)
      self._nuIdxr = self._nu:genIndexer()

      self.volTermFunc = function(w, dx, idx, q, out)
         return ConstDiffusion["volTermVarInSpace"](self, w, dx, idx, q, out)
      end
      self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return ConstDiffusion["surfTermVarInSpace"](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
      self.boundarySurfTermFunc = function(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
         return ConstDiffusion["boundarySurfTermVarInSpace"](self, dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
      end
   else
      self.volTermFunc = function(w, dx, idx, q, out)
         return ConstDiffusion["volTermConstInSpace"](self, w, dx, idx, q, out)
      end
      self.surfTermFunc = function(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
         return ConstDiffusion["surfTermConstInSpace"](self, dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
      end
      self.boundarySurfTermFunc = function(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
         return ConstDiffusion["boundarySurfTermConstInSpace"](self, dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
      end
   end
end

-- Methods.
function ConstDiffusion:numEquations() return 1 end
function ConstDiffusion:numWaves() return 1 end
function ConstDiffusion:isPositive(q)
   if q[1] > 0.0 then
      return true
   else
      return false
   end
end

-- Flux in direction dir.
function ConstDiffusion:flux(dir, qIn, fOut)
   assert(false, "ConstDiffusion:flux: NYI!")
end

-- Riemann problem for ConstDiffusion equation.
function ConstDiffusion:rp(dir, delta, ql, qr, waves, s)
   assert(false, "ConstDiffusion:rp: NYI!")
end

-- Compute q-fluctuations.
function ConstDiffusion:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "ConstDiffusion:qFluctuations: NYI!")
end

-- Maximum wave speed.
function ConstDiffusion:maxSpeed(dir, w, dx, q)
   assert(false, "ConstDiffusion:maxSpeed: NYI!")
   return 0.0
end

-- Volume integral term for use in DG scheme.
function ConstDiffusion:volTerm(w, dx, idx, q, out)
   return self.volTermFunc(w, dx, idx, q, out)
end
function ConstDiffusion:volTermConstInSpace(w, dx, idx, q, out)
   local cflFreq = self._volTerm(w:data(), dx:data(), self._nu:data(), q:data(), out:data())
   return cflFreq
end
function ConstDiffusion:volTermVarInSpace(w, dx, idx, q, out)
   self._nu:fill(self._nuIdxr(idx), self._nuPtrl)
   local cflFreq = self._volTerm(w:data(), dx:data(), self._nuPtrl:data(), q:data(), out:data())
   return cflFreq
end

-- Surface integral term for use in DG scheme.
function ConstDiffusion:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   return self.surfTermFunc(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
end
function ConstDiffusion:surfTermConstInSpace(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end
function ConstDiffusion:surfTermVarInSpace(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self._nu:fill(self._nuIdxr(idxl), self._nuPtrl)
   self._nu:fill(self._nuIdxr(idxr), self._nuPtrr)
   self._surfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), self._nuPtrl:data(), self._nuPtrr:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end

-- Contribution from surface integral term at the boundaries for use in DG scheme.
function ConstDiffusion:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
   return self.boundarySurfTermFunc(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
end
function ConstDiffusion:boundarySurfTermConstInSpace(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
   self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), edge, self._nu:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end
function ConstDiffusion:boundarySurfTermVarInSpace(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, edge, ql, qr, outl, outr)
   self._nu:fill(self._nuIdxr(idxl), self._nuPtrl)
   self._nu:fill(self._nuIdxr(idxr), self._nuPtrr)
   self._boundarySurfTerms[dir](wl:data(), wr:data(), dxl:data(), dxr:data(), idxl:data(), idxr:data(), edge, self._nuPtrl:data(), self._nuPtrr:data(), ql:data(), qr:data(), outl:data(), outr:data())
   return 0
end

function ConstDiffusion:setAuxFields(auxFields)
end


return ConstDiffusion
