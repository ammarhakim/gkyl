-- Gkyl ------------------------------------------------------------------------
--
-- vlasov equation using Hamiltonian formulation
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CanonicalModDecl = require "Eq.pbData.CanonicalModDecl"
local Proto = require "Lib.Proto"
local HamiltonianBase = require "Eq.HamiltonianBase"
local ProjectOnBasis = require "Updater.ProjectOnBasis"

-- start from HamiltonianBase class
local HamilVlasov = Proto(HamiltonianBase)

-- ctor
function HamilVlasov:init(tbl)
   -- initialize generic Hamiltonian equation
   HamilVlasov.super.init(self, tbl)
   assert(self._confBasis, "HamilVlasov: must specify confBasis")

   local charge = assert(tbl.charge, "HamilVlasov: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "HamilVlasov: must specify mass using 'mass' ")
   self._qbym = charge/mass -- only q/m ratio is ever needed

   self._hasPhi = tbl.hasPhi
   self._isElectromagnetic = assert(tbl.hasA==false, "HamilVlasov: electromagnetic not yet implemented!")

   self._ndim = self._basis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._ndim - self._cdim

   -- store pointers to C kernels implementing volume and surface terms
   local nm, p = self._basis:id(), self._basis:polyOrder()
   self._volTerm = CanonicalModDecl.selectVol(nm, self._cdim, self._vdim, p)
   self._surfTerms = CanonicalModDecl.selectSurf(nm, self._cdim, self._vdim, p)
   self._disContCorrectionSurfTerms = CanonicalModDecl.selectDisContCorrectionSurf(nm, self._cdim, self._vdim, p)

   self:initHamilTimeIndep()  -- equation dependent
end

function HamilVlasov:initHamilTimeIndep()
   local projectKE = ProjectOnBasis {
      onGrid = self._grid,
      basis = self._basis,
      evaluate = function (t, xn)
                    local ke = 0
                    for vi=self._cdim+1, self._ndim do
                       ke = ke + xn[vi]*xn[vi]/2
                    end
                    return ke
                 end
                    
   }
   projectKE:advance(0, {}, {self.hamilTimeIndep})
end

function HamilVlasov:setAuxFields(auxFields)
   if self._hasPhi then 
      -- get phi, and scale by q/m
      self.phi = auxFields[1].phi
      self.phi:scale(self._qbym)

      -- time-dependent part of hamiltonian is q/m*phi
      self:setHamiltonian(self.phi)
   end
end

return HamilVlasov
