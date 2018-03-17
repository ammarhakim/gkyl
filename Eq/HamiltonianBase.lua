-- Gkyl ------------------------------------------------------------------------
--
-- Generic Hamiltonian equation base class
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local ffi = require "ffi"
local xsys = require "xsys"

local Hamiltonian = Proto()

-- ctor
function Hamiltonian:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.grid, "Hamiltonian equation must specify a grid")
   self._basis = assert(tbl.basis, "Hamiltonian equation must specify a basis")

   -- initialize hamiltonian
   self.hamiltonian = DataStruct.Field {
     onGrid = self._grid,
     numComponents = self._basis:numBasis(),
     ghost = {1,1}
   }
   self.hamiltonian:clear(0.0)

   -- initialize time-dependent and time-independent parts of hamiltonian
   self.hamilTimeDep = self.hamiltonian
   self.hamilTimeIndep = self.hamiltonian

   -- make pointers and indexers
   self.hamPtr = self.hamiltonian:get(1)
   self.hamPtrL, self.hamPtrR = self.hamPtr, self.hamPtr
   self.hamIdxr = self.hamiltonian:genIndexer()
end

function Hamiltonian:setHamiltonian(hamilTimeDepIn)
   -- check if hamilTimeDep is on same grid as self.hamiltonian
   -- if not, project hamilTimeDep onto grid of self.hamiltonian 
   local hamilTimeDep
   if hamilTimeDepIn:compatible(self.hamiltonian) then
      hamilTimeDep = hamilTimeIndepIn
   else
      assert(false, "Using hamiltonian terms on different grids not yet implemented!")
   end

   -- set hamiltonian
   self.hamiltonian:combine(1.0, self.hamilTimeIndep, 1.0, hamilTimeDep)
end

-- Volume integral term for use in DG scheme
function Hamiltonian:volTerm(w, dx, idx, q, out)
   self.hamiltonian:fill(self.hamIdxr(idx), self.hamPtr)
   return self._volTerm(w:data(), dx:data(), q:data(), self.hamPtr:data(), out:data())
end
-- Surface integral term for use in DG scheme
function Hamiltonian:surfTerm(dir, w, dx, maxs, idxl, idxr, ql, qr, outl, outr)
   self.hamiltonian:fill(self.hamIdxr(idxl), self.hamPtrL)
   self.hamiltonian:fill(self.hamIdxr(idxr), self.hamPtrR)
   return self._surfTerms[dir](w:data(), dx:data(), ql:data(), qr:data(), self.hamPtrL:data(), self.hamPtrR:data(), outl:data(), outr:data())
end

return Hamiltonian
