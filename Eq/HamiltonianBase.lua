-- Gkyl ------------------------------------------------------------------------
--
-- Generic Hamiltonian equation base class
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local DataStruct = require "DataStruct"
local ConfToPhase = require "Updater.ConfToPhase"
local EqBase = require "Eq.EqBase"
local xsys = require "xsys"

local Hamiltonian = Proto(EqBase)

-- ctor
function Hamiltonian:init(tbl)
   -- get grid and basis
   self._grid = assert(tbl.onGrid, "Hamiltonian equation must specify a grid")
   self._basis = assert(tbl.basis or tbl.phaseBasis, "Hamiltonian equation must specify a basis")
   self._confBasis = tbl.confBasis
   self._positivity = xsys.pickBool(tbl.positivity, false)

   self._ndim = self._grid:ndim()

   -- set hamiltonian discontinuity direction flags
   self._hamilDisCont = {}
   for d = 1, self._ndim do
      self._hamilDisCont[d] = false
   end
   if tbl.hamilDisContDirs then
      for i, d in ipairs(tbl.hamilDisContDirs) do
         self._hamilDisCont[d] = true
      end
   end

   -- initialize hamiltonian
   self.hamiltonian = DataStruct.Field {
     onGrid = self._grid,
     numComponents = self._basis:numBasis(),
     ghost = {1,1}
   }
   self.hamiltonian:clear(0.0)

   -- set up time-independent part of hamiltonian
   self.hamilTimeIndep = DataStruct.Field {
     onGrid = self._grid,
     numComponents = self._basis:numBasis(),
     ghost = {1,1}
   }
   self.hamilTimeIndep:clear(0.0)

   self.hamilTimeDep = DataStruct.Field {
     onGrid = self._grid,
     numComponents = self._basis:numBasis(),
     ghost = {1,1}
   }
   self.hamilTimeDep:clear(0.0)

   -- make pointers and indexers
   self.hamPtr = self.hamiltonian:get(1)
   self.hamPtrL, self.hamPtrR = self.hamiltonian:get(1), self.hamiltonian:get(1)
   self.hamIdxr = self.hamiltonian:genIndexer()

   -- updater to project a configuration space field to phase space
   if self._confBasis then
      self.accumulateConfToPhase = ConfToPhase {
         onGrid = self._grid,
         confBasis = self._confBasis,
         phaseBasis = self._basis,
         operation = "accumulate"
      }
      self.assignConfToPhase = ConfToPhase {
         onGrid = self._grid,
         confBasis = self._confBasis,
         phaseBasis = self._basis,
         operation = "assign"
      }
   end
 
   self.ioFrame = 0
end

function Hamiltonian:setHamiltonian(hamilTimeDepIn)
   -- set hamiltonian by summing time-dependent and time-independent parts
   -- check if hamilTimeDep is on same grid as self.hamiltonian
   -- if not, project hamilTimeDep onto grid of self.hamiltonian
   if hamilTimeDepIn:compatible(self.hamiltonian) then
      -- set hamiltonian = hamilTimeIndep + hamilTimeDep
      self.hamiltonian:combine(1.0, self.hamilTimeIndep, 1.0, hamilTimeDepIn)
   else
      -- set hamiltonian = hamilTimeIndep
      self.hamiltonian:combine(1.0, self.hamilTimeIndep)
      -- accumulate hamiltonian += hamilTimeDep
      self.accumulateConfToPhase:advance(0, 0, {1.0, hamilTimeDepIn}, {self.hamiltonian})
      self.assignConfToPhase:advance(0, 0, {1.0, hamilTimeDepIn}, {self.hamilTimeDep})
   end
end

-- Volume integral term for use in DG scheme
function Hamiltonian:volTerm(w, dx, idx, q, out)
   self.hamiltonian:fill(self.hamIdxr(idx), self.hamPtr)
   return self._volTerm(w:data(), dx:data(), self.hamPtr:data(), q:data(), out:data())
end
-- Surface integral term for use in DG scheme
function Hamiltonian:surfTerm(dir, cfl, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr)
   self.hamiltonian:fill(self.hamIdxr(idxl), self.hamPtrL)
   self.hamiltonian:fill(self.hamIdxr(idxr), self.hamPtrR)
   if self._hamilDisCont[dir] then 
      self._disContCorrectionSurfTerms[dir](wr:data(), dxr:data(), maxs, self.hamPtrL:data(), self.hamPtrR:data(), ql:data(), qr:data(), outl:data(), outr:data()) 
   end
   return self._surfTerms[dir](cfl, wr:data(), dxr:data(), maxs, self.hamPtrR:data(), ql:data(), qr:data(), outl:data(), outr:data())
end

function Hamiltonian:writeHamiltonian(distIo, tm)
   distIo:write(self.hamiltonian, string.format("%s_%d.bp", "hamil", self.ioFrame), tm)
   self.ioFrame = self.ioFrame + 1
end

return Hamiltonian
